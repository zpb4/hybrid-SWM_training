#setwd('z:/oro_nonstat/')
library(ranger)
library(lime)

print(paste('start',Sys.time()))

#specifications
hym_site<-'ORO'
vers<-'err13'
noise_reg<-T  #use noise regularized coefficients?
seed<-10

#load predictor arrays
hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('./analysis_data/hym_predmat_4c_',hym_site,'.rds',sep=''))
rf_err_corr<-readRDS(paste('./model_output/hymod_rf-err-corr_cal_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

ix<-seq(as.Date('1988-10-01'),as.Date('2018-09-30'),'day')
ix2<-as.POSIXlt(ix)

idx_cal<-which(ix=='1988-10-01'):which(ix=='1998-09-30')
idx_val<-which(ix=='1998-10-01'):which(ix=='2004-09-30')
idx_trn<-which(ix=='1988-10-01'):which(ix=='2004-09-30')
idx_tst<-which(ix=='2004-10-01'):which(ix=='2018-09-30')

idx_cal<-which(ix=='1988-10-01'):which(ix=='1989-09-30')

rf_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0','err -1','err -2','err -3')))

pred_df<-as.data.frame(hym_predmat_hist[idx_tst,rf_idx])
pred_df_4c<-as.data.frame(hym_predmat_4c[idx_tst,rf_idx])

#lime setup
explainer <- lime(pred_df,rf_err_corr)
explainer_4c <- lime(pred_df_4c,rf_err_corr)

explanation <- explain(pred_df, explainer,n_labels=1,n_features = length(rf_idx))
explanation_4c <- explain(pred_df_4c, explainer_4c,n_labels=1,n_features = length(rf_idx))

saveRDS(explanation,paste('./model_output/explanation_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
saveRDS(explanation_4c,paste('./model_output/explanation_4c_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))

feat_wt<-matrix(explanation$feature_weight,ncol=length(rf_idx),byrow=T)
feat_wt_4c<-matrix(explanation_4c$feature_weight,ncol=length(rf_idx),byrow=T)

lime_feat_wt<-feat_wt
colnames(lime_feat_wt)<-colnames(hym_predmat_hist[,rf_idx])
rownames(lime_feat_wt)<-rownames(hym_predmat_hist[idx_tst,])

lime_feat_wt_4c<-feat_wt_4c
colnames(lime_feat_wt_4c)<-colnames(hym_predmat_4c[,rf_idx])
rownames(lime_feat_wt_4c)<-rownames(hym_predmat_4c[idx_tst,])

saveRDS(lime_feat_wt,paste('./model_output/lime_feat_wt_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))
saveRDS(lime_feat_wt_4c,paste('./model_output/lime_feat_wt_4c_hymod_v-',vers,'_nreg=',noise_reg,'.rds',sep=''))

print(paste('end',Sys.time()))

rm(list=ls());gc()


##########################################END################################