# hybrid-SWM training and familiarization
Overview: This file has some redundancy with the associated WaterProgramming blog post [here](https://waterprogramming.wordpress.com/2024/03/11/nonstationary-stochastic-watershed-modeling/). It is intended to
describe elements of the hybrid-SWM implementation in more detail than a surface-level implementation in the README file. The individual scripts are also heavily commented to describe the workings of each element of the model.   
### Overview of the GRRIEN repository structure
In the Steinschneider group, we cover elements of the GRRIEN repository setup as part of internal training [here](https://github.com/SteinschneiderLab/lab-manual/tree/main/training/open_research). See also Rohini's
blog post on the same [here](https://waterprogramming.wordpress.com/2023/03/06/introducing-the-grrien-analysis-framework-defining-standards-for-reproducible-and-robust-supervised-learning-of-earth-surface-processes-at-large-spatial-scales/), which describe Elizabeth Carter's paper (Syracuse University, former Steinschneider Group member) that formalizes this method.
![image info](figures_tables/GRRIEn.png "GRRIEN")
#### _GRRIEN repository_   
   
The README file and structure of this repository has been constructed to follow closely with this GRRIEN structure, excepting some elements that are more appropriate for large geospatial data problems. As noted in the group training on this subject, this is one of many possible formulations of a good repository structure and useful as a starting point to build out your own preferred approach.
### Data processing
The data processing for the hybrid SWM (data_process.R) is relatively straightforward. This script simply tidies up the state-variable and simulation timeseries files from the raw_data folder and calculates lag-1 to 3 errors for fitting the Random Forest (RF) error correction model. The data generation via the hydrologic models (SAC-SMA and HYMOD) for the Feather River at Oroville (ORO) site used for this example is not covered in detail in this training. You can consult the manuscript and supporting information files associated with this training for more details on the actual setup and optimization of the hydrologic models, which was all done by Sungwook Wi. The specifics of this are not super important to the main idea. What is important is the setup of the 'model-as-truth' experimental design used to generate these data. See the associated blog post for a high-level description of this approach. The figures showing the study watershed and describing the experimental design are shown below for convenience:
![image info](figures_tables/fig1.png "Study area")   
#### _Study area_      
   
![image info](figures_tables/fig2.png "Model-as-truth experimental design")
#### _Model-as-truth experimental design_   
### Model fitting
Fitting of the hybrid-SWM is done via a staged, two-step process. The first step is an error correction model that is fit as a predictive model between the state-variables and the raw predictive errors, including lag-1 to 3 errors to account for autocorrelation. The second step is a dynamic residual model (DRM) that is a time-varying, state-variable dependent distributional model for the residuals of the error correction step. We'll go through the steps in more detail below with reference to the 'model_train.R' script. Here is the figure from the manuscript describing the model fitting procedure for reference:
![image info](figures_tables/fig3.png "hybrid SWM")
#### _Hybrid-SWM_  

### Random Forest error correction model
As mentioned in the model fitting section, the error correction model is simply a predictive model between the state-variables and the raw errors. We implement this step using a Random Forest (RF) model, due to their robustness to overfitting. RF models are a relatively early ML technique that leverage an original methodology called Classification and Regression Trees (CART). Fitting of individual CART trees IS prone to overfitting, but RF builds and ensemble of trees using a few randomization techniques (bagging, random feature selection). The output of the RF is the 'majority vote' of the ensemble of trees and coupled with the other randomization techniques, prevents overfitting to the training data. See Rohini's post [here](https://waterprogramming.wordpress.com/2019/02/04/intro-to-machine-learning-part-5-bagging/) for further details.
![image info](figures_tables/RF.png "hybrid SWM")
#### _Random Forest_  

Importantly, the RF error correction model is fitted to a calibration subset of the training data, leaving the validation subset for fitting of the DRM. To interpret variable importance, we can use the output of the fitted RF model directly, which calculated feature importance aggregated over all the trees in the ensemble. This output is described in more detail in the 'RF variable importance' figure and text below.   
   
We can also employ LIME (Local, Interpretable, Model-agnostic Explanation) as a form of explainable AI (xAI) to the fitted RF model to understand local feature importance, down to the granularity of individual timesteps. [Here](docs/LIME.pdf) is a short pictorial depiction of LIME. We describe the LIME procedure as implemented in this work ('lime' package in R) in the associated section below.

### Dynamic residual model (DRM)
The DRM is a little bit methodologically dense. It embeds linear models for all the parameters of a distribution, the skew exponential power (SEP) distribution, inside a log-likelihood function that can be maximized against the empirical residuals in the validation subset after RF error correction. This allows the DRM to capture 'out-of-sample' residual uncertainty and time-varying, state-variable dependent properties of those residuals. We employ a noise regularization technique to the MLE procedure that helps smooth out the estimation of the coefficients against the data. As noted, the SEP distribution is the heart of this model and most importantly, the parameterizations for kurtosis ($\beta$) and skewness ($\xi$) shown . What this means practically is that the result of the DRM fitting procedure ('model_train.R') yields a set of coefficients (9 state-variable coefficients + 1 intercept) for each SEP parameter's linear model. So there is a state-variable dependent linear model for the stdev term $\sigma_t$, a linear model for the kurtosis parameter $\beta_t$, etc. Importantly, there is no model for the mean $\mu$ because the residuals are assumed to be debiased by the error correction procedure. Thus, given any timeseries of state-variables, the fitted linear models will generate a timeseries of the 4 parameters of the DRM ($\sigma_t, \beta_t, \xi_t, \phi_t$). Functionally, this means that generation of new residual at time = $t$ is a random sample from the $SEP(0,\sigma_t,\beta_t,\xi_t)$ modified by the residual at $t-1$  ($\epsilon_{t-1}$) through the $\phi_t$ lag-1 autocorrelation term.  
![image info](figures_tables/sep.png "SEP distribution")
#### _SEP distribution_  

## Code implementation
This is all quite a lot to parse. Let's take a look at the code implementation, which may clarify some of these higher level descriptions. The model scripts included in this repo have been commented quite extensively, so I'll just try to highlight key sections of the code to explain what they are doing.   
### data_process
As noted above, the data processing part of the procedure is not super hard to follow. The script is just downloading the raw data, applying labels to it, and doing some simple offsetting to establish lag 0 to lag 3 versions of all state-variables and the error timeseries. The state variables used for HYMOD are shown in the table below. The matrix output by the data-processing script is repeated from lag 0 to lag 3, but the only lagged information used in the actual modeling is the lag 1 to 3 errors. Note how the errors are defined in the script as the 'truth' model minus the 'process' model outputs (L66-71).
![image info](figures_tables/statvar_table.png "State variables")
#### _State variables_ 
### model_train
This script fits both the RF error correction model and the DRM. Up until L35, the script sets up the date/time indices for the calibration, validation, and test periods and organizes the predictor matrices for both the RF error correction (all lag 0 state variables + lag 1:3 errors) and the DRM (all lag 0 state variables. The actual fitting of the RF error correction model requires only one line of code on L42 from the 'ranger' package. There is no hyperparameter selection for the RF model as we found it not to be super helpful, so the RF model uses the default settings of 500 trees and $\sqrt{K}$ for the random feature selection parameter 'mtry'. This model is fit on the calibration subset of the training data and then in L48-56, the RF model is used to debias the errors of the validation set for fitting of the DRM.   

Fitting of the DRM model takes up a bit more of the script than the RF model and is targeted to the RF debiased residuals $\epsilon_t$ in the validation subset. L65-94 define the bounded constraints for the optimization problem. While $\beta_t$ and $\log_{10}(\xi_t)$ can take values relatively symmetric about zero, the $\sigma_t$ and $\phi_t$ should be positive and no greater than 1 for $\phi_t$. Most notably, we found that constraining all linear coefficients for $\sigma_t$ to be > 0  and also setting the minimum intercept to the lowest decile of the $|\epsilon_t|$ ensured stability and achieved good MLE solutions. L96-112 deal with normalization of the predictors, where predictors are scaled and centered for $\beta_t, \xi_t$, but scaled and referenced to a minimum value of zero for $\sigma_t, \phi_t$. These normalization vectors are saved for use later in SWM generataion. This specific setup was somewhat undesirable from an ease of implementation standpoint, but enabled more stable and consistent MLE results. Finally, a technique called 'noise regularization' is used in L114-132 to help smooth out the MLE estimation. This procedures adds random Gaussian noise to both the predictors and the predictand distribution. We specifically add $\mathcal{N}(0,0.2)$ to the predictors and $\mathcal{N}(0,0.1)$ to the distribution of $\epsilon_t$ per recommendations of Rothsfuss et al. (2019).   
   
MLE is accomplished with a gradient based optimization routine from the R 'optimx' package where all 40 parameters (9 state variable coefficients and 1 intercept term per parameter) are fit simultaneously. L151-172 configure the coefficient matrices for later use and print the log-likelihood values for diagnosis of the results. The optimization can be sensitive and in the actual manuscript implementation, we use multiple seeded runs and select the best seed (by log-likelihood value) to use in generation. This is not implemented here for simplicity, but you can change the seed or optimization technique if needed. If the optimization does not find solution, it often spits out the starting parameter values and will generally achieve log-likelihood values in the negative thousands (e.g. -6000). Good solutions are typically in the positive range. The log-likelihood function is configured to return extremely large negative values if internal function constraints are not met (e.g. NA values, Inf values, etc).
   
![image info](figures_tables/drm_functions.PNG "DRM functions")
#### _DRM functions_ 
### model_simulate
This script simulates from the fitted model and is combination of 3 subscripts, 'model_simulate_test', 'model_simulate_4c', and 'model_simulate_benchmark'. The '_test' and '_4c' scripts are functionally identical except that one uses the historical forcings and the other uses the 4C warmed forcings. The '_benchmark' script is a 'static' SWM model that is described in more detail in the manuscript (and will not be discussed in detaill here). If you actually want to run the scripts, I recommend running each separately, particularly if you have multiple compute nodes at your disposal. The simulate script is not fast (~30-45 min) on HPC for 100X samples since it has an autoregressive term in the RF error correction model and must simulate each timestep in sequence. This requires an instantiation of the RF model in 'predict' mode at each timestep, which is time consuming. Also, the scripts are configured to generate across the entired timespan of the data even thought the Test/Test+4C subset is what is of most interest for the study results. There are pre-staged .rds files run from these exact scripts that you can use in lieu of generating them yourselves.   

L8-16 setup the parallel implementation of the simulation using the 'foreach' loop structure. The script is configured to use all available compute cores (with 2 left free) to run the simulation. L18-65 are basic setup information for the simulation, while L71-74 apply the normalization to the new data using the scaling and centering vectors from the 'model_train' script. It is important that the normalization be implemented this way and not 'refitted' to the new data, because that could introduce unintended biases to the model predictions. L79-114 is the core error generation code, and specifically L99-111 is where the core simulation occurs. The script describes it pretty well in the comments. The simulation model first generates a set of residuals from the fitted DRM which are then integrated into the simulation via the RF error correction model in 'predict' mode. Importantly, the RF error correction model uses the previous three timesteps of errors to predict the next timestep, so these errors must be sequentially loaded into the prediction matrix at each timestep. The math for this simulation procedure is below, where the tilde indicates a simulated value and $\theta_{SV,t}$ is the vector of state-variables at time = $t$:   
   
$\tilde{e_t}=f(\theta_{SV,t},\tilde{e_{t-3}},\tilde{e_{t-2}},\tilde{e_{t-1}})+\tilde{\epsilon_t}$   

After generation of the errors, the SWM is simulated by adding these errors back to the original process model simulation and ensuring that any negative flow simulations are corrected to zero (L120-124). The error generation procedure is repeated as many times as the desired number of simulation samples (100 by default in these scripts).

### model_lime
This script simulates

![image info](figures_tables/fig5.png "Error correction result")
#### _Error correction residuals_  

### RF variable importance figure (Figure 6)
This is the variable importance data output directly from the RF implementation used in this work ('ranger' package in R). This 
![image info](figures_tables/fig6.png "RF variable importance")
#### _RF variable importance_  

Below is the result of the LIME procedure applied to subsets of the data to accentuate features that contribute most to inferring changing biases between the Test and Test+4C scenarios.   
![image info](figures_tables/fig7.png "LIME")
#### _LIME_  

The application of the DRM is shown for the same periods, demonstrating its ability to adapt in the Test and Test+4C case.  
![image info](figures_tables/fig8.png "DRM fits")
#### _Fitted DRM result_  
   
Finally, interpretability of the model can be relatively easily accomplished looking at properties of the parameters, such as their seasonality.
![image info](figures_tables/fig9.png "DRM seasonality")
#### _DRM parameter seasonality_  

Comparison of the hybrid SWM to a 'static' approach noting the inability of the static model to emulate changes in the Test+4C period.
![image info](figures_tables/fig10.png "hybrid SWM vs static SWM")
#### _Hybrid SWM vs static SWM_  

Timeseries plot in 2011 showing implications of shifts in biases between Test and Test+4C case. Of note during snowmelt periods and flow recession.
![image info](figures_tables/fig11.png "hybrid SWM timeseries plot")
#### _Hybrid SWM timeseries plot_  

Comparison of hybrid SWM and static SWM in terms of emulating extreme high and low flows in Test and Test+4C.
![image info](figures_tables/fig12.png "hybrid SWM vs static SWM - high/low flow extremes")
#### _Hybrid SWM vs static SWM - high/low flow extremes_ 
