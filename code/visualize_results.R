#visualization of nonstationary SWM
#setwd('z:/hybrid-SWM_training/')

#Figure scripts below are reference to figures in the manuscript


#Plot Figure 4 - Test and Test+4C comparison of raw errors for '5wet' and '5dry' periods
#brief demo of how you can set-up individual scripts as plotting functions
source('./code/functions/plot_functions/plot_fig4.R')
x = 14
plot_fig4(14)

rm(list=ls());gc()

#Note: a better and more generalizable way to do this is to define basic plotting functions
#that you can arrange how you like in the script and save the output.

#For instance, for figure 4 plot here, would be better to have a 'boxplot_error_comparison' script
#that you could call with different data entries for the Test or Test+4C scenarios and then maybe
#a 'plot_hydrograph' kind of script for panel c

#For brevity, the remainder of the plotting scripts are just sourced:

#Plot Figure 5 - Test and Test+4C comparison of raw errors and RF error-corrected residuals
source('./code/functions/plot_functions/plot_fig5.R')

rm(list=ls());gc()

#Plot Figure 6 - Variable importance for the RF model
source('./code/functions/plot_functions/plot_fig6.R')

rm(list=ls());gc()

#Plot Figure 7 - LIME analysis of the RF model for March subset (can change month in the script)
source('./code/functions/plot_functions/plot_fig7.R')

rm(list=ls());gc()

#Plot Figure 8 - Residual analysis for 3 months comparing empirical residuals to fitted DRM
source('./code/functions/plot_functions/plot_fig8.R')

rm(list=ls());gc()

#Plot Figure 9 - DRM parameter timeseries plots
source('./code/functions/plot_functions/plot_fig9.R')

rm(list=ls());gc()

#Plot Figure 10 - Error comparison boxplots between hybrid SWM and static SWM for Test and Test+4C
source('./code/functions/plot_functions/plot_fig10.R')

rm(list=ls());gc()

#Plot Figure 11 - hybrid SWM timeseries plots for 6 month period in 2011
source('./code/functions/plot_functions/plot_fig11.R')

rm(list=ls());gc()

#Plot Figure 12 - High and Low flow quantile comparison between hybrid SWM and static SWM for Test and Test+4C
source('./code/functions/plot_functions/plot_fig12.R')

rm(list=ls());gc()

#Plot Table 2 - DRM coefficient table
#NOTE: Could not verify this was working on my machine due to some issues with the R packages due to updates
source('./code/functions/plot_functions/plot_table2.R')

rm(list=ls());gc()



###########################################END###############################################