#setwd('z:/oro_nonstat/')
library(fGarch)
library(ranger)
library(doParallel)

print(paste('start',Sys.time()))

#parallelization setup
parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

#model specifications
n<-100 #number of samples to generate
ar<-3
hym_site<-'ORO'
samp_type<-'split' # 'split' 'skip'
gen_period<-'hist-all' # 'test' 'cal' 'val' 'train' 'hist-all'
vers<-'err13'
noise_reg<-T  #use noise regularized coefficients?
seed<-10

#load data
source('./code/functions/GL_maineqs.R')

#DRM fitted coefficients
dyn_res_coef<-readRDS(paste('./model_output/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
#zero ref for predictor variables from training
pred_mat_zero<-readRDS(paste('./model_output/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
#historical predictor array
hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))
#normalization vector from fitted period
norm_vec<-readRDS(paste('./model_output/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))
#fitted RF error correction model
rf_err_corr<-readRDS(paste('./model_output/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#array label indices
rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))
                                                                       
#do not include lag error terms in dynamic residual prediction
res_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))


#define date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

if(samp_type=='split'){
  idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
  idx_hist<-which(ix=='1988-10-01'):which(ix=='2018-09-30')
  idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
  idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
  idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')
}

if(gen_period=='train'){gen_idx<-idx_trn}
if(gen_period=='cal'){gen_idx<-idx_cal}
if(gen_period=='val'){gen_idx<-idx_val}
if(gen_period=='test'){gen_idx<-idx_tst}
if(gen_period=='hist-all'){gen_idx<-idx_hist}

#predictor matrix for the RF error correction model
pred_mat<-hym_predmat_hist[,rf_idx]
pred_mat_rf<-pred_mat[gen_idx,]

#define and normalize predictor matrix for the DRM (no lag terms)
dyn_res_preds<-hym_predmat_hist[,res_idx]
pred_mat_scale<-(dyn_res_preds[gen_idx,]-matrix(rep(norm_vec[1,],length(gen_idx)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(gen_idx)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(gen_idx)),ncol=dim(pred_mat_scale)[2],byrow=T)

#empirical error vector 
e_t<-hym_predmat_hist[gen_idx,'err 0']

#model simulation code written in parallel foreach loop
syn_out<-foreach(m = 1:n,.combine='cbind',.packages=c('fGarch','ranger'),.inorder=F) %dopar% {
  #generate random residuals from the DRM
  syn_res<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=e_t)
  #output vector
  syn_gen<-rep(0,(length(gen_idx)))
  
  #lengthen output vector by number of ar terms and add random samples from residual vector as first '1:ar' inputs
  if(ar>0){
    syn_gen<-rep(0,(length(gen_idx)+ar));syn_gen[1:ar]<-sample(syn_res,ar)
  }
  #copy RF predictor matrix
  dat<-pred_mat_rf
  #define indices for the lagged error terms
  if(ar>0){
    ar_terms<-paste('err',0:-ar)
    rf_gen_idx<-sort(which((colnames(pred_mat_rf)%in%ar_terms)==T))
    }
  
  #main generation code that sequences through length of data
  for(i in 1:length(gen_idx)){
    #with non-zero ar terms, select errors from lag-1 to lag-'ar' and put in predictor matrix for RF
    if(ar>0){
      #idx is lag-1 to lag-'ar'
      ar_idx<-(i+ar-1):i
      #put lag terms in predictor matrix for timestep i
      dat[i,rf_gen_idx]<-syn_gen[ar_idx] #ar3
    }
    #first predict the error for timestep i
    err<-predict(rf_err_corr,data=dat[i,])$predictions
    #then add that prediction to the residual for timestep i; put in output matrix at i+ar (offset removed at end)
    syn_gen[i+ar]<-err+syn_res[i]
  }
  #return the output vector minus the first 1:ar terms that were used to start the generation
  return(syn_gen[(ar+1):length(syn_gen)])
}

#save the generated errors
saveRDS(syn_out,paste('./model_output/hymod_syn-err_',hym_site,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'_samps=',n,'.rds',sep=''))

#to generate SWM output, need to add errors to the hydro model simulation
sim<-matrix(rep(hym_predmat_hist[gen_idx,'sim 0'],n),ncol=n,byrow=F)
swm_out<-sim+syn_out

#ensure any negative flow simulations are corrected to zero
swm_out[swm_out<0]<-0

#save SWM data
saveRDS(swm_out,paste('./model_output/hymod_syn-flow_',hym_site,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'_samps=',n,'.rds',sep=''))
print(paste('end',Sys.time()))

rm(list=ls());gc()

##############################END##################