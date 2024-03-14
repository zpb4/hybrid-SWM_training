#visualization of nonstationary SWM
#setwd('z:/hybrid-SWM_training/')

#Figure scripts below are reference to figures in the manuscript


#Plot Figure 4 - Test and Test+4C comparison of raw errors for '5wet' and '5dry' periods
#brief demo of how you can set-up individual scripts as plotting functions
source('./code/functions/plot_functions/plot_fig4.R')
x = 14
plot_fig4(14)

#Note: a better and more generalizable way to do this is to define basic plotting functions
#that you can arrange how you like in the script and save the output.

#For instance, for figure 4 plot here, would be better to have a 'boxplot_error_comparison' script
#that you could call with different data entries for the Test or Test+4C scenarios and then maybe
#a 'plot_hydrograph' kind of script for panel c

#For brevity, the remainder of the plotting scripts are just sourced:

#Plot Figure 5 - Test and Test+4C comparison of raw errors and RF error-corrected residuals
source('./code/functions/plot_functions/plot_fig5.R')