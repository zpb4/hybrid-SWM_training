#setwd('z:/oro_nonstat/')
library(ranger)
library(fGarch)
source('./code/functions/mm-cfs_conversion.R')

hym_site<-'ORO'
sma_site<-'ORO'
vers<-'err13'
noise_reg<-T  #use noise regularized coefficients?
seed=10

#wy subset to compare
#'5wet': c(2005,2006,2011,2016,2017)
#'5dry': c(2007,2008,2012,2014,2015)
#'tst-all' : 2005:2018
wy_comp<-c(2005,2006,2011,2016,2017)
wy_tag<-'5wet'

hym_kcfs_conv<-as.numeric(area_calc[hym_site,2])*mm_to_cfs/1000
sma_kcfs_conv<-as.numeric(area_calc[sma_site,2])*mm_to_cfs/1000

sma_hist_vars<-readRDS(paste('./analysis_data/sma_hist_vars_',sma_site,'.rds',sep=''))
hym_hist_vars<-readRDS(paste('./analysis_data/hym_hist_vars_',sma_site,'.rds',sep=''))
hym_4c_vars<-readRDS(paste('./analysis_data/hym_4c_vars_',sma_site,'.rds',sep=''))
hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('./analysis_data/hym_predmat_4c_',hym_site,'.rds',sep=''))
sma_predmat_hist<-readRDS(paste('./analysis_data/sma_predmat_hist_',sma_site,'.rds',sep=''))
sma_predmat_4c<-readRDS(paste('./analysis_data/sma_predmat_4c_',sma_site,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')

ix_tst<-seq(as.Date('2004-10-01'),as.Date('2018-09-30'),'day')
ixx_tst<-as.POSIXlt(ix_tst)
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

obs_tst<-sma_hist_vars[idx_tst,'obs']
tst_wy<-2005:2018
obs_sum<-c()

for(i in 1:length(tst_wy)){
  obs_sum[i]<-sum(obs_tst[which(ix_tst==as.Date(paste(tst_wy[i]-1,'-10-01',sep=''))):which(ix_tst==as.Date(paste(tst_wy[i],'-09-30',sep='')))])
}

srt_obs<-sort(obs_sum,index.return=T)
tst_wy[tail(srt_obs$ix,5)] #5 wettest
#[1] 2016 2005 2011 2006 2017
tst_wy[head(srt_obs$ix,5)] #5 driest
#[1] 2014 2015 2008 2007 2012

ix_comp<-c()
idx_comp<-c()
for(i in 1:length(wy_comp)){
  ix_comp<-c(ix_comp,seq(as.Date(paste(wy_comp[i]-1,'-10-01',sep='')),as.Date(paste(wy_comp[i],'-09-30',sep='')),'day'))
  idx_comp<-c(idx_comp,which(ix==paste(wy_comp[i]-1,'-10-01',sep='')):which(ix==paste(wy_comp[i],'-09-30',sep='')))
}
ix_comp<-as.Date(ix_comp,origin = '1970-01-01')
ixx_comp<-as.POSIXlt(ix_comp)

rf_err_corr<-readRDS(paste('./model_output/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
norm_vec<-readRDS(paste('./model_output/hymod_norm-vec_hist_',hym_site,'.rds',sep=''))


dyn_res_coef<-readRDS(paste('./model_output/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
pred_mat_zero<-readRDS(paste('./model_output/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

#--------------------------------------------------------------
rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))

#do not include lag error terms in dynamic residual prediction
res_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))
#plot raw errors vs corrected errors
db_err<-predict(rf_err_corr,data=hym_predmat_hist[idx_comp,rf_idx])$predictions
err_db<-hym_predmat_hist[idx_comp,'err 0']-db_err

db_err_4c<-predict(rf_err_corr,data=hym_predmat_4c[idx_comp,rf_idx])$predictions
err_db_4c<-hym_predmat_4c[idx_comp,'err 0']-db_err_4c
#--------------------------------------------------------------------------------------
#plot AR residuals vs predicted residuals
source('./code/functions/GL_maineqs.R')

#test
pred_mat<-hym_predmat_hist[,rf_idx]

dyn_res_preds<-hym_predmat_hist[,res_idx]
pred_mat_rf<-pred_mat[idx_comp,]
pred_mat_scale<-(dyn_res_preds[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds)[2]),byrow=T)
pred_mat_zero_min<-pred_mat_scale-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale)[2],byrow=T)

set.seed(seed)
syn_res<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min,beta_var=pred_mat_scale,xi_var=pred_mat_scale,phi_var=pred_mat_zero_min,et=err_db)*hym_kcfs_conv
syn_res[is.na(syn_res)==T]<-0


#test4c
pred_mat_4c<-hym_predmat_4c[,rf_idx]

dyn_res_preds_4c<-hym_predmat_4c[,res_idx]
pred_mat_rf_4c<-pred_mat_4c[idx_comp,]
pred_mat_scale_4c<-(dyn_res_preds_4c[idx_comp,]-matrix(rep(norm_vec[1,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T))/matrix(rep(norm_vec[2,],length(idx_comp)),ncol=(dim(dyn_res_preds_4c)[2]),byrow=T)
pred_mat_zero_min_4c<-pred_mat_scale_4c-matrix(rep(pred_mat_zero,length(idx_comp)),ncol=dim(pred_mat_scale_4c)[2],byrow=T)

set.seed(seed)
syn_res_4c<-syn_gen_mv_ar1_lin(dyn_res_coef,sig_var=pred_mat_zero_min_4c,beta_var=pred_mat_scale_4c,xi_var=pred_mat_scale_4c,phi_var=pred_mat_zero_min_4c,et=err_db_4c)*hym_kcfs_conv
syn_res_4c[is.na(syn_res_4c)==T]<-0

brks<-rbind(c(-100,seq(-5,5,0.4),100),c(-100,seq(-5,5,0.4),100),c(-100,seq(-5,5,0.4),100))
xlm<-rbind(c(-4,4),c(-4,4),c(-4,4))
x_axis<-rbind(seq(-4,4,2),seq(-4,4,2),seq(-4,4,2))
yupr<-c(1.25,1.25,1.25)
pad<-c(0.015,0.015,0.015)
mths<-c(2,3,4)

mth<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

png(paste('./figures_tables/fig8.png',sep=''),width=768,height=384)
par(mfrow=c(2,3),oma=c(0,5,0,0),mar=c(1,4,3,0.5),mgp=c(2.5,0.5,0),tcl=-0.2,cex.lab=2,cex.axis=2,cex.main=3,font.main=4)

for(i in 1:3){
  seas<-which(ixx_comp$mon==(mths[i]-1))
  kd_syn_err<-density(syn_res[seas],n=512)
  if(i==1){dens_lab='Density'}else{dens_lab=''}
  hist(err_db[seas]*hym_kcfs_conv,breaks=brks[i,],xlim=xlm[i,],ylim=c(0,yupr[i]),freq=F,main=mth[mths[i]],col='skyblue',xlab='',ylab=dens_lab,axes=F)
  axis(1,at=x_axis[i,],labels=F)
  if(i>=1){axis(2,at=seq(0,3,.5),labels=seq(0,3,.5))}
  #if(i>1){axis(2,at=seq(0,3,.5),labels=F)}
  kd_y<-kd_syn_err$y
  kd_y[kd_y<=pad[i]]<-pad[i]
  lines(kd_syn_err$x,kd_y,col='red',lwd=3)
  if(i==1){mtext('Test',side=2,cex=2.5,line=5)
    legend('topleft',c('emp','sim'),col=c('skyblue','red'),lwd=c(6,3),cex=2,bty='n')}
}

par(mar=c(4,4,1,0.5))
    
for(i in 1:3){
  seas<-which(ixx_comp$mon==(mths[i]-1))
  kd_syn_err<-density(syn_res_4c[seas],n=512)
  if(i==1){dens_lab='Density'}else{dens_lab=''}
  hist(err_db_4c[seas]*hym_kcfs_conv,breaks=brks[i,],xlim=xlm[i,],ylim=c(0,yupr[i]),freq=F,main='',col='goldenrod1',xlab='Residual (kcfs)',ylab=dens_lab,axes=F)
  axis(1,at=x_axis[i,],labels=x_axis[i,])
  if(i>=1){axis(2,at=seq(0,3,.5),labels=seq(0,3,0.5))}
  #if(i>1){axis(2,at=seq(0,3,.5),labels=F)}
  kd_y<-kd_syn_err$y
  kd_y[kd_y<=pad[i]]<-pad[i]
  lines(kd_syn_err$x,kd_y,col='red',lwd=3)
  if(i==1){mtext('Test+4C',side=2,cex=2.5,line=5)
    legend('topleft',c('emp','sim'),col=c('goldenrod1','red'),lwd=c(6,3),cex=2,bty='n')}
}

dev.off()

#################################################END#########################################
