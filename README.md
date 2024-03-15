# hybrid-SWM_training
Practical training exercise repository for hybrid-SWM implementation from Brodeur et al. (2024) 'A hybrid, non-stationary Stochastic Watershed Model (SWM) for uncertain hydrologic simulations under climate change'. Supports Waterprogramming (Reed Group) blog post [here](https://waterprogramming.wordpress.com/2024/03/11/nonstationary-stochastic-watershed-modeling/). This README describes the basic components of the repository and how to implement the code. It also and serves as an example of a README file that includes the recommended components of a GRRIEN repository. The 'training_exercise.md' covers the algorithmic implementation of the model in detail for training and familiarization purposes.
### Repository overview   
This repository contains most of the standard elements of a GRRIEN repository, excluding the setup.sh functionality, 'data_source', and 'temp_data' elements that are used for more advanced repository structures setup to source large geospatial datasets:
- README: this file
- training_exercise.md: the file that describes in detail the algorithmic implementation of the hybrid SWM
- run_code.cmd: file to run code elements individually on an HPC with SLURM using the 'sbatch' command, i.e. 'sbatch run_code.cmd'
- raw_data: repository for the raw state-variable and simulation information for the SAC-SMA 'truth' model and the HYMOD 'process' model as described in the blog post. Data are included for the Feather River at Oroville (ORO) watershed.
- analysis_data: repository to store pre-processed data outputs from the 'data_process.R' code
- model_output: repository to store model objects and outputs from the 'model_train.R' and 'model_simulate.R' codefiles
    - NOTE: 100x sample runs and LIME results are archived here from previous generation runs; simulate code is lengthy to run and you can use these pre-staged runs for the visualizations
- figures_tables: repository to store figures generated from the 'visualize_results.R' code
- code: repository for all codefiles
    - functions: subrepository of 'code' that stores basic functions and routines for the primary codefiles
- docs: documentation associated with this repo, primarily the manuscript and supporting info
### Code implementation instructions
#### 1. Run ./code/setup.R   
   _Description: Installs required R packages_
#### 2. Run ./code/data_process.R   
   _Description: Preprocesses the raw data from .txt files in R data structure arrays_       
   I: raw data .txt files   
   O: preprocessed data arrays     
      
   _RUNTIME: <1 min_ 
#### 3. Run ./code/model_train.R
   _Description: Fits both components of the hybrid SWM to the training data_      
   I: pre-processed state-variable and hydrologic model simulation data from the 'data_process' step   
   O: trained model objects saved as .rds files to the 'model_output' folder   
      
   _RUNTIME: <5 min_
#### 4. Run ./code/model_simulate.R   
   _Description: Simulates from the fitted hybrid SWM and a 'static' benchmark SWM model for both the historical and 4C warmed cases across the entire range of the data (WY1989-2018). Must specify no. of samples to simulate (default is 100)._   
   I: trained model objects and pre-processed data arrays   
   O: SWM simulations saved as .rds arrays   
      
   _RUNTIME: >1 hour on HPC for 100x samples; if just wanting to familiarize with model and code, you can use the pre-staged simulation data in the 'model_outputs' folder._  

   _NOTE: This script combines 'model_simulate_test.R', 'model_simulate_test4c.R', and 'model_simulate_benchmark.R' for illustrative purposes. Recommend running these scripts individually if you'd like to speed up generation, particularly on HPC resources where each script can be assigned to a different node. The 'model_simulate_benchmark.R' runs in <5 min._   
#### 5. Run ./code/model_lime.R   
   _Description: Runs the Local Interpretable Model Agnostic (LIME) routine against the trained Random Forest error correction model_   
   I: trained RF model object and pre-processed data arrays   
   O: LIME explanation object      
      
   _RUNTIME: <30 min; if just wanting to familiarize with model and code, you can use the pre-staged LIME data in the 'model_outputs' folder._   
#### 6. Run ./code/visualize_results.R
   _Description: Generates a set of figures described in more detail in the training_exercise.md file_   
   I: trained model objects, pre-processed data arrays, SWM simulations   
   O: figures/tables saved as .png files      

### Dataset citation
Brodeur, Z., Wi, S., Shabestanipour, G., Lamontagne, J. R., & Steinschneider, S. A hybrid, non-stationary Stochastic Watershed Model (SWM) for uncertain hydrologic simulations under climate change. Water Resources Research, in review.
### Contact information
Zach Brodeur (Cornell University), zpb4@cornell.edu
### Preferred citation
Brodeur, Z. P. (2024). hybrid-SWM_training: Mar 11, 2024 release (v1.0.0) [Software]. Github. https://github.com/zpb4/hybrid-SWM_training
### Change request policy
Please raise a GitHub issue on this site or contact the repository owner via email
### License
GNU general public license 3.0. See 'LICENSE.md' file for details
