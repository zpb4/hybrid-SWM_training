_Note: Repo is currently under construction; should be complete by 3/12/24_
# hybrid-SWM training and familiarization
Overview: This file has some redundancy with the associated WaterProgramming blog post [here](https://waterprogramming.wordpress.com/2024/03/11/nonstationary-stochastic-watershed-modeling/). It is intended to
describe elements of the hybrid-SWM implementation in more detail than a surface-level implementation in the README file. 
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
As mentioned in the model fitting section, the error correction model is simply a predictive model between the state-variables and the raw errors. We implement this step using a Random Forest (RF) model, due to their robustness to overfitting. RF models are a relatively early ML technique that leverage an original methodology called Classification and Regression Trees (CART). Fitting of individual CART trees IS prone to overfitting, but RF builds and ensemble of trees using a few randomization techniques (bagging, random feature selection). The output of the RF is the 'majority vote' of the ensemble of trees and coupled with the other randomization techniques, prevents overfitting to the training data.
![image info](figures_tables/RF.png "hybrid SWM")
#### _Random Forest_  

The code implementation of the RF model is straightforward, as is prediction from the fitted model. Importantly, the RF error correction model is fitted to a calibration subset of the training data, leaving the validation subset for fitting of the DRM. To interpret variable importance, we can use the output of the fitted RF model directly, which calculated feature importance aggregated over all the trees in the ensemble. 
![image info](figures_tables/fig6.png "RF variable importance")
#### _RF variable importance_  

We can also employ LIME (Local, Interpretable, Model-agnostic Explanation) as a form of explainable AI (xAI) to the fitted RF model to understand local feature importance, down to the granularity of individual timesteps. [Here](docs/LIME.pdf) is a short pictorial depiction of LIME. Below is the result of the LIME procedure applied to subsets of the data to accentuate features that contribute most to inferring changing biases between the Test and Test+4C scenarios.   
![image info](figures_tables/fig7.png "LIME")
#### _LIME_  

### Random Forest error correction model
