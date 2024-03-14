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

#//////////////////////////////////////////////////////////////////////////////////////////////
#Test+4C generation
#NOTE: all elements are the same as above, but use the 4C warmed inputs and simulations
#comments are 

print(paste('start 4C',Sys.time()))

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

dyn_res_coef<-readRDS(paste('./model_output/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
pred_mat_zero<-readRDS(paste('./model_output/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#-------------------------------------------------------------------------------

pred_mat<-hym_predmat_4c[,rf_idx]
dyn_res_preds<-hym_predmat_4c[,res_idx]

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

saveRDS(syn_out,paste('./model_output/hymod_syn-err_',hym_site,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'_samps=',n,'.rds',sep=''))

sim<-matrix(rep(hym_predmat_4c[gen_idx,'sim 0'],n),ncol=n,byrow=F)
swm_out<-sim+syn_out
swm_out[swm_out<0]<-0
saveRDS(swm_out,paste('./model_output/hymod_syn-flow_',hym_site,'_',gen_period,'_',samp_type,'_nreg=',noise_reg,'_',n,'X_v-',vers,'_samps=',n,'.rds',sep=''))
print(paste('end 4C',Sys.time()))


#////////////////////////////////////////////////////////////////////////////////////////////////
#3. Fit and simulate from benchmark model
#NOTE: this model is a simpler implementation following Schoups & Vrugt (2010) that does not have the hybrid and dynamic capabilities of the primary model
#It is used as a comparison in some of the generated figures

library(fGarch)

#specifications
seed<-10
hym_site<-'ORO'

#load predictor arrays
hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('./analysis_data/hym_predmat_4c_',hym_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_hist<-which(ix=='1988-10-01'):which(ix=='2018-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

ix_trn<-seq(as.Date('1988-10-01'),as.Date('2004-09-30'),'day')
ixx_trn<-as.POSIXlt(ix_trn)

ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)

ix_hist<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ixx_hist<-as.POSIXlt(ix_hist)

#fit benchmark model for error correction
err_hist<-hym_predmat_hist[idx_trn,'err 0']
err_hist_db<-err_hist

err_mns<-c()

for(i in 1:12){
  seas<-which(ixx_trn$mon==(i-1))
  err_mns[i]<-mean(err_hist[seas])
  err_hist_db[seas]<-err_hist[seas]-rep(err_mns[i],length(seas))
}

saveRDS(err_mns,'./model_output/hymod_benchmark_err-mns.rds')

#GL-SGED model
sim_trn<-hym_predmat_hist[idx_trn,'sim 0']
sim_hist<-hym_predmat_hist[,'sim 0']
sim_4c<-hym_predmat_4c[,'sim 0']

source('./code/functions/GL_maineqs_uv.R')

lb<-c(0.0001,0,0,0,0,-0.99,0.1) #lower bounds for sig0, sig1, phi1, phi2, phi3, beta, xi; -1 for beta gives an error
ub<-c(20,5,1,1,1,5,10) #upper bounds for sig0, sig1, phi1, phi2, phi3, beta, xi
st<-c(.5,.5,.5,.5,.5,0,1) #starting parameters for sig0, sig1, phi1, phi2, phi3, beta, xi

gl_par<-array(NA,c(12,7))

for(i in 1:12){
  seas<-which(ixx_trn$mon==(i-1))
  et_sort<-sort(abs(err_hist_db[seas]))
  min_sig<-max(0.0001,mean(et_sort[1:round(0.1*length(seas))]))
  lb[1]<-min_sig
  gl_mle<-optim(par=st,fn=GL_fun_noscale_ar3_bfgs,inflow=sim_hist[seas],et=err_hist_db[seas],
                method = 'L-BFGS-B',lower = lb, upper = ub,
                control = list(fnscale=-1,maxit=100000))
  gl_par[i,]<-gl_mle$par
}

saveRDS(gl_par,'./model_output/hymod_benchmark_gl-par.rds')

#hist
syn_err_hist<-array(NA,c(length(ix_hist),n))
syn_flow_hist<-array(NA,c(length(ix_hist),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_hist$mon==(i-1))
    et<-et_syn(sim_hist[seas],sim_hist[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err<-et+rep(err_mns[i],length(seas))
    syn_flow<-sim_hist[seas]+syn_err
    syn_flow[syn_flow<0]<-0
    syn_err<-syn_flow - sim_hist[seas]
    syn_err_hist[seas,m]<-syn_err
    syn_flow_hist[seas,m]<-syn_flow
  }
}

saveRDS(syn_err_hist,'./model_output/hymod_benchmark_syn-err_hist.rds')
saveRDS(syn_flow_hist,'./model_output/hymod_benchmark_syn-flow_hist.rds')

#4c
syn_err_4c<-array(NA,c(length(ix_hist),n))
syn_flow_4c<-array(NA,c(length(ix_hist),n))

for(m in 1:n){
  for(i in 1:12){
    seas<-which(ixx_hist$mon==(i-1))
    et<-et_syn(sim_4c[seas],sim_4c[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    syn_err<-et+rep(err_mns[i],length(seas))
    syn_flow<-sim_4c[seas]+syn_err
    syn_flow[syn_flow<0]<-0
    syn_err<-syn_flow - sim_4c[seas]
    syn_err_4c[seas,m]<-syn_err
    syn_flow_4c[seas,m]<-syn_flow
  }
}

saveRDS(syn_err_4c,'./model_output/hymod_benchmark_syn-err_4c.rds')
saveRDS(syn_flow_4c,'./model_output/hymod_benchmark_syn-flow_4c.rds')

rm(list=ls());gc()

##############################END##################