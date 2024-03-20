#fit model between hymod and sacsma
#setwd('d:/hybrid-SWM_training/')
library(fGarch)
library(ranger)
library(optimx)

#specifications
vers<-'err13' 
seed<-10
hym_site<-'ORO'
sma_site<-'ORO'

#load predictor arrays
hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))
sma_hist_vars<-readRDS(paste('./analysis_data/sma_hist_vars_',sma_site,'.rds',sep=''))

rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))

#do not include lag error terms in dynamic residual prediction
res_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

#----------------------------------------------------
#1) Prepare data
#raw errors are target for RF error correction model
rf_err<-hym_predmat_hist[,'err 0']
#pull 'rf_idx' indices from larger prediction matrix
rf_pred<-hym_predmat_hist[,rf_idx]

#---------------------------------------------------
#2) Fit Random Forest error correction model
#RF model is fit under default settings (no hyperparameter tuning)
#For the purposes of error modeling, we did not find much utility in hyperparameter tuning (see manuscript for details)
set.seed(seed)
rf_err_corr<-ranger(x=rf_pred[idx_cal,],y=rf_err[idx_cal],importance = 'impurity')
saveRDS(rf_err_corr,paste('./model_output/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#----------------------------------------------------
#3) Fit multivariate, dynamic residual model (DRM)

#first apply the RF error correction model in the validation period and debias the errors
#this yields residuals (epsilon_t) that the DRM is fit to

#3a.debias calibration errors with RF error correction model

#predict from RF error correction model
db_val<-predict(rf_err_corr,data=rf_pred[idx_val,])$predictions
#debiased errors are simply raw error minus the RF error correction model predictions
err_db_val<-rf_err[idx_val]-db_val

#3b. Setup constrained optimization problem for the DRM

#source the log-likelihood function
#see GL_maineqs and GL_subeqs code for more details
source('./code/functions/GL_maineqs.R')

#define lower bound for sigma intercept
srt_err<-sort(abs(err_db_val)) #use absolute value of errors as a proxy for sigma_t to set upper and lower bounds
min_sig<-mean(srt_err[1:round(length(srt_err)*0.1)])
max_sig<-max(abs(err_db_val))

#define bounds for sigma [par1], beta [par2], xi [par3], and phi [par4]
#starting parameter values
st_arr<-array(NA,c(4,(length(res_idx)+1)))
st_arr[1,]<-c(0.5,rep(0.5,length(res_idx)))
st_arr[2,]<-c(0,rep(0,length(res_idx)))
st_arr[3,]<-c(0,rep(0,length(res_idx)))
st_arr[4,]<-c(0.5,rep(0,length(res_idx)))

#lower bound arrays
lb_arr<-array(NA,c(4,(length(res_idx)+1)))
lb_arr[1,]<-c(min_sig,rep(0,length(res_idx))) #sigma intercept bounded at a 'min_sig' value, must be >0
lb_arr[2,]<-c(-0.99,rep(-5,length(res_idx)))
lb_arr[3,]<-c(-1,rep(-1,length(res_idx)))
lb_arr[4,]<-c(0,rep(-1,length(res_idx)))

#upper bound arrays
ub_arr<-array(NA,c(4,(length(res_idx)+1)))
ub_arr[1,]<-c(max_sig,rep(max_sig,length(res_idx))) #constrain upper sigma bounds to maximum 
ub_arr[2,]<-c(5,rep(5,length(res_idx)))
ub_arr[3,]<-c(1,rep(1,length(res_idx)))
ub_arr[4,]<-c(1,rep(1,length(res_idx)))

#arrays above must be vectorized (flattened) to go into optimization setup
start<-as.vector(st_arr)
lb<-as.vector(lb_arr)
ub<-as.vector(ub_arr)

#3c.Normalize predictors
#select columns for DRM (state variables, no lag terms)
pred_mat<-hym_predmat_hist[,res_idx]
#center and scale the predictors
pred_mat_scale<-scale(pred_mat[idx_trn,])
ctr_vec<-attributes(pred_mat_scale)['scaled:center']$'scaled:center'
scale_vec<-attributes(pred_mat_scale)['scaled:scale']$'scaled:scale'

#3d.Retain normalization vectors to use in simulation 
#save normalization vector to reuse on Test and Test+4C data
norm_vec<-rbind(ctr_vec,scale_vec)
saveRDS(norm_vec,paste('./model_output/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))

#Test normalization vector on fitted data for later use; should return 0 if matrices are the same
#check with manual calculation
pred_mat_comp<-(pred_mat[idx_trn,]-matrix(rep(ctr_vec,length(idx_trn)),ncol=dim(pred_mat)[2],byrow=T))/matrix(rep(scale_vec,length(idx_trn)),ncol=dim(pred_mat)[2],byrow=T)
print(max(abs(pred_mat_comp-pred_mat_scale))) #should return 0

#3e. Noise regularization

#Noise regularization is a relatively straightforward technique to add Gaussian noise to predictors (normalized state variables) and 
#predictand (residual distribution) to get a smoother and more generalizable MLE estimation of parameters

#Noise regularization description:
#inputs/predictand -> z* = z + z*norm(0,sd) [Klotz et al. 2022, p16]
#recommended setting [Rothsfuss et al. 2019, p15] sd_x=0.2, sd_y=0.1 
#So, add white noise of N(0,0.2) to predictors and N(0,0.1) to residual distribution

#Add noise to predictor matrix
noise_pred_mat<-pred_mat_scale + pred_mat_scale * matrix(rnorm(length(idx_trn)*length(res_idx),sd=0.2),nrow=length(idx_trn))
#Add noise to residuals
err_db_val<-err_db_val + err_db_val * rnorm(length(err_db_val),sd=0.1)

#create a zero reference predictor matrix for the sigma and phi coefficients (both should be >0), save zero reference vector
noise_pred_mat_zero<-apply(noise_pred_mat,2,min)
saveRDS(noise_pred_mat_zero,paste('./model_output/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
noise_pred_mat_zero_min<-noise_pred_mat-matrix(rep(noise_pred_mat_zero,length(idx_trn)),ncol=dim(noise_pred_mat)[2],byrow=T)

#final prediction matrices
dyn_res_preds<-noise_pred_mat[idx_val,]
dyn_res_preds_zero<-noise_pred_mat_zero_min[idx_val,]

#3f. Optimize the DRM via MLE
#NOTE: MLE for this many parameters is sensitive and sometimes does not find solutions. In manuscript implementation, we use
#multiple seeded runs and choose the best solution by the log-likelihood. This takes some time though. See notes below if issues getting a solutions
noise_dyn_res_mle<-optimx(par=start,fn=GL_fun_mv_ar1_lin,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,
                          xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val,noise=rep(0,length(err_db_val)),neg=F,
                          lower = lb,upper = ub,method = 'nlminb',
                          control = list(maximize=T,all.methods=F))

#NOTE: If issues with getting an MLE solution (i.e. the 'fitted' coefficients are the start coefficients above),
#try different seeds or different optimization method below:
#method options: ('lbfgsb','nlminb','spg','Rcgmin','Rvmmin','bobyqa','nmkb','hjkb')
#can also run 'all.methods=T' in control list to try out all 8 options in a single run (takes a while)

#retrieve the fitted DRM coefficients from the optimization object
noise_dyn_res_coef<-c()
for(i in 1:length(start)){
  noise_dyn_res_coef[i]<-noise_dyn_res_mle[[i]][1]
}
noise_dyn_res_coef

#input the parameters into a coefficient matrix
noise_coef_mat<-matrix(noise_dyn_res_coef,nrow=4,byrow=F)
colnames(noise_coef_mat)<-c('intcpt',res_idx)
rownames(noise_coef_mat)<-c('sigma','beta','xi','phi')

#run the log-likelihood function to verify the log-likelihood (ll) values
#Note: the baseline 'start' parameters will achieve a ll somewhere in the negative thousands (e.g -6000)
#Good solutions should be greater than 0 and are often in the positive thousands (e.g. 3000)
noise_ll<-GL_fun_mv_ar1_lin(noise_dyn_res_coef,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val,noise=rep(0,length(err_db_val)),neg=F)
#calculate parameter timeseries for each of the 4 DRM parameters
noise_param_out<-GL_fun_mv_ar1_lin_params(noise_dyn_res_coef,sig_var=dyn_res_preds_zero,beta_var=dyn_res_preds,xi_var=dyn_res_preds,phi_var=dyn_res_preds_zero,et=err_db_val)

#print the log-likelihood values and coefficient matrix for real-time analysis (helpful for HPC implementation)
print(noise_ll)
print(noise_coef_mat)

#save coefficient matrices and parameter timeseries
saveRDS(noise_dyn_res_coef,paste('./model_output/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
saveRDS(noise_coef_mat,paste('./model_output/hymod_noise-dyn-res-coef-mat_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
saveRDS(noise_param_out,paste('./model_output/hymod_noise-dyn-res-params_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

rm(list=ls());gc()

##########################################END#################################################
