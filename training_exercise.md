_Note: Repo is currently under construction; should be complete by 3/12/24_
# hybrid-SWM training and familiarization
Overview: This file has some redundancy with the associated WaterProgramming blog post [here](https://waterprogramming.wordpress.com/2024/03/11/nonstationary-stochastic-watershed-modeling/). It is intended to
describe elements of the hybrid-SWM implementation in more detail than surface-level implementation in the README file. 
### Overview of the GRRIEN repository structure
In the Steinschneider group, we cover elements of the GRRIEN repository setup as part of internal training [here](https://github.com/SteinschneiderLab/lab-manual/tree/main/training/open_research). See also Rohini's
blog post on the same [here](https://waterprogramming.wordpress.com/2023/03/06/introducing-the-grrien-analysis-framework-defining-standards-for-reproducible-and-robust-supervised-learning-of-earth-surface-processes-at-large-spatial-scales/), which describe Elizabeth Carter's paper (Syracuse University, former Steinschneider Group member) that formalizes this method.
![image info](figures_tables/GRRIEn.png "GRRIEN")
The README file and structure of this repository has been constructed to follow closely with this GRRIEN structure, excepting some elements that are more appropriate for large geospatial data problems. As noted in the group training on this subject, this is one of many possible formulations of a good repository structure and useful as a starting point to build out your own preferred approach.
### Data processing
The data processing for the hybrid SWM (data_process.R) is relatively straightforward. This script simply tidies up the state-variable and simulation timeseries files from the raw_data folder and calculates lag-1 to 3 errors for fitting the Random Forest (RF) error correction model. The data generation via the hydrologic models (SAC-SMA and HYMOD) for the Feather River at Oroville (ORO) site used for this example is not covered in detail in this training. You can consult the manuscript and supporting information files associated with this training for more details on the actual setup and optimization of the hydrologic models, which was all done by Sungwook Wi. The specifics of this are not super important to the main idea. What is important is the setup of the 'model-as-truth' experimental design used to generate these data. See the associated blog post for a high-level description of this approach. The figures showing the study watershed and describing the experimental design are shown below for convenience:
![image info](figures_tables/fig1.png "Study area")
![image info](figures_tables/fig2.png "Model-as-truth experimental design")
### Model fitting
