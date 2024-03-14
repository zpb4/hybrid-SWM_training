#setwd('z:/oro_nonstat/')
library(gt)
library(stringr)

hym_site<-'ORO'
sma_site<-'ORO'
#specifications
vers<-'err13' 
noise_reg=T
seed=10

hym_predmat_hist<-readRDS(paste('./analysis_data/hym_predmat_hist_',hym_site,'.rds',sep=''))

res_idx<-sort(which(colnames(hym_predmat_hist)%in%c('precip 0','tavg 0','sim 0','runoff 0','baseflow 0','et 0','swe 0','upr_sm 0','lwr_sm 0')))

statvar_names<-colnames(hym_predmat_hist[,res_idx])
statvar_names<-str_remove(statvar_names,' 0')
#plot AR residuals vs predicted residuals

dyn_res_coef<-readRDS(paste('./model_output/hymod_noise-dyn-res-coef_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))
pred_mat_zero<-readRDS(paste('./model_output/hymod_noise-pred-mat-zero_val_',hym_site,'_v-',vers,'_seed',seed,'.rds',sep=''))

tab <- matrix(as.numeric(dyn_res_coef), ncol=4, byrow=T)
row_rearrange<-c('intcpt','sim','runoff','baseflow','upr_sm','lwr_sm','swe','et','tavg','precip')
row_def <- c('intcpt',statvar_names[1:10])
inp_idx<-c()
for(i in 1:length(row_rearrange)){
  inp_idx[i]<-which(row_def==row_rearrange[i])
}
tab_inp<-tab[inp_idx,]

tab_mat<-cbind(row_rearrange,tab_inp)
colnames(tab_mat) <- c('statvar','sigma','beta','xi','phi')
#rownames(tab_mat) <- row_rearrange
#tab_out<-data.frame(tab_mat,row.names = row_rearrange)
tab_out<-as.data.frame(tab_mat)
tab_out[,2:5]<-tab_inp

tab_gt<-gt(tab_out,rowname_col = 'statvar')%>%
  tab_header(
    title = "Residual Model Coefficients",
  )%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = statvar
  ))%>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(
      columns = c(statvar,sigma,beta,xi,phi)
    )
  )%>%
  fmt_number(
    columns=c(sigma,beta,xi,phi),
    decimals = 3
  )%>%
  cols_label(
    sigma = html("&sigma;<sub>t</sub>"),
    beta = html("&beta;<sub>t</sub>"),
    xi = html("log<sub>10</sub>(&xi;<sub>t</sub>)"),
    phi = html("&phi;<sub>t</sub>"),
    )%>%
    cols_align(
        align = 'center',
    )

gtsave(tab_gt,paste('./figures_tables/table2.png',sep=''))


