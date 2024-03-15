#'Static' SWM model
#setwd('z:/oro_nonstat/')
library(fGarch)

#specifications
seed<-10
hym_site<-'ORO'

#load predictor arrays
hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))
hym_predmat_4c<-readRDS(paste('./analysis_data/hym_predmat_4c_',hym_site,'.rds',sep=''))

#define date/time indices
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

#retrieve raw errors
err_hist<-hym_predmat_hist[idx_trn,'err 0']
err_hist_db<-err_hist

#-----------------------------------------------
#Implement GL-SEP model for static SWM per manuscript
#1) define monthly means to debiase errors
err_mns<-c()

for(i in 1:12){
  seas<-which(ixx_trn$mon==(i-1))
  err_mns[i]<-mean(err_hist[seas])
  err_hist_db[seas]<-err_hist[seas]-rep(err_mns[i],length(seas))
}

saveRDS(err_mns,'./model_output/hymod_benchmark_err-mns.rds')

#-------------------------------------------------------------
#2) Fit static SWM model to errors
#GL-SEP model error model from Schoups & Vrugt (2010)
#Maximize a generalized likelihood function with linear heteroscedastic transform and AR(3) model embedded in likelihood function

#require simulations for heteroscedastic modeling
sim_trn<-hym_predmat_hist[idx_trn,'sim 0']
sim_hist<-hym_predmat_hist[,'sim 0']
sim_4c<-hym_predmat_4c[,'sim 0']

#retrieve univariate (uv) equations
source('./code/functions/GL_maineqs_uv.R')

#constrained optimization initial settings
lb<-c(0.0001,0,0,0,0,-0.99,0.1) #lower bounds for sig0, sig1, phi1, phi2, phi3, beta, xi; -1 for beta gives an error
ub<-c(20,5,1,1,1,5,10) #upper bounds for sig0, sig1, phi1, phi2, phi3, beta, xi
st<-c(.5,.5,.5,.5,.5,0,1) #starting parameters for sig0, sig1, phi1, phi2, phi3, beta, xi

gl_par<-array(NA,c(12,7))

#optimize GL-SEP error model by month
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

#---------------------------------------------
#3) Generate from fitted model for Test (historical) scenario
#This script actually generates across the entire length of the data
n<-100

syn_err_hist<-array(NA,c(length(ix_hist),n))
syn_flow_hist<-array(NA,c(length(ix_hist),n))

#generation script
for(m in 1:n){
  #generate by month from the monthly fitted parameter sets
  for(i in 1:12){
    seas<-which(ixx_hist$mon==(i-1))
    #simulate residuals
    et<-et_syn(sim_hist[seas],sim_hist[seas],gl_par[i,1],gl_par[i,2],gl_par[i,3],gl_par[i,4],gl_par[i,5],gl_par[i,6],gl_par[i,7])[[1]]
    #add in the estimated mean error bias
    syn_err<-et+rep(err_mns[i],length(seas))
    #swm output is errors + simulations
    syn_flow<-sim_hist[seas]+syn_err
    #correct simulations then modify errors to ensure flows >= 0
    syn_flow[syn_flow<0]<-0
    syn_err<-syn_flow - sim_hist[seas]
    syn_err_hist[seas,m]<-syn_err
    syn_flow_hist[seas,m]<-syn_flow
  }
}

saveRDS(syn_err_hist,'./model_output/hymod_benchmark_syn-err_hist.rds')
saveRDS(syn_flow_hist,'./model_output/hymod_benchmark_syn-flow_hist.rds')

#---------------------------------------------
#4) Generate from fitted model for Test+4C (historical) scenario
#As above, this script actually generates across the entire length of the data
n<-100

syn_err_4c<-array(NA,c(length(ix_hist),n))
syn_flow_4c<-array(NA,c(length(ix_hist),n))

#generation script as above except with Test+4C data
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

##########################################END#################################################

