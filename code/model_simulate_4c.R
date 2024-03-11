#setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(doParallel)

print(paste('start',Sys.time()))

parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

n<-10
ar<-3
hym_site<-'ORO'
samp_type<-'split' # 'split' 'skip'
gen_period<-'4c-all' # 'test4c' 'cal4c' 'val4c' 'train4c' '4c-all'
vers<-'err13'
noise_reg<-T #use noise regularized coefficients?

hym_predmat_4c<-readRDS(paste('./analysis_data/hym_predmat_4c_',hym_site,'.rds',sep=''))
norm_vec<-readRDS(paste('./model_output/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))
rf_err_corr<-readRDS(paste('./model_output/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))


rf_idx<-sort(which(colnames(hym_predmat_4c)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))
#do not include lag error terms in dynamic residual prediction
res_idx<-sort(which(colnames(hym_predmat_4c)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

if(samp_type=='split'){
  idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
  idx_hist<-which(ix=='1988-10-01'):which(ix=='2018-09-30')
  idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
  idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
  idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')
}

if(gen_period=='train4c'){gen_idx<-idx_trn}
if(gen_period=='cal4c'){gen_idx<-idx_cal}
if(gen_period=='val4c'){gen_idx<-idx_val}
if(gen_period=='test4c'){gen_idx<-idx_tst}
if(gen_period=='4c-all'){gen_idx<-idx_hist}

#load data
source('./code/functions/GL_maineqs.R')
if(noise_reg==T){
  dyn_res_coef<-readRDS(paste('./model_output/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('./model_output/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
}
if(noise_reg==F){
  dyn_res_coef<-readRDS(paste('./model_output/hymod_dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
  pred_mat_zero<-readRDS(paste('./model_output/hymod_pred-mat-zero_val_',hym_site,'.rds',sep=''))
}

#-------------------------------------------------------------------------------

pred_mat<-hym_predmat_4c[,rf_idx]
dyn_res_preds<-hym_predmat_4c[,res_idx]

#generate
pred_mat_rf<-pred_mat[gen_idx,]
pred_mat_scale<-(dyn_res_preds[gen_idx,]-matrix(rep(norm_vec[1,],length(gen_idx)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(gen_idx)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(gen_idx)),ncol=dim(pred_mat_scale)[2],byrow=T)

e_t<-hym_predmat_4c[gen_idx,'err 0']

syn_out<-foreach(m = 1:n,.combine='cbind',.packages=c('fGarch','ranger'),.inorder=F) %dopar% {
  
  syn_res<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=e_t)
  syn_gen<-rep(0,(length(gen_idx)))
  if(ar>0){
    syn_gen<-rep(0,(length(gen_idx)+ar));syn_gen[1:ar]<-sample(syn_res,ar)
  }
  
  dat<-pred_mat_rf
  if(ar>0){
    ar_terms<-paste('err',0:-ar)
    rf_gen_idx<-sort(which((colnames(pred_mat_rf)%in%ar_terms)==T))
    }
  
  for(i in 1:length(gen_idx)){
    if(ar>0){
      ar_idx<-(i+ar-1):i
      dat[i,rf_gen_idx]<-syn_gen[ar_idx] #ar3
    }
    err<-predict(rf_err_corr,data=dat[i,])$predictions
    syn_gen[i+ar]<-err+syn_res[i]
  }
  
  return(syn_gen[(ar+1):length(syn_gen)])
}

saveRDS(syn_out,paste('./model_output/hymod_syn-err_',hym_site,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))

sim<-matrix(rep(hym_predmat_4c[gen_idx,'sim 0'],n),ncol=n,byrow=F)
swm_out<-sim+syn_out
swm_out[swm_out<0]<-0
saveRDS(swm_out,paste('./model_output/hymod_syn-flow_',hym_site,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'.rds',sep=''))
print(paste('end',Sys.time()))

rm(list=ls());gc()

##############################END##################