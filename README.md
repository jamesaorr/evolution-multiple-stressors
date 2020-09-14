# evolution-multiple-stressors
Data and code used for all statistical analyses and figures for the paper: "Rapid evolution generates synergism between multiple stressors and complicates ecological restoration". 

Code for analysis and figure/video construction are in three R notebooks, which should be run in order: 

- **1_Evolution-Multiple-Stressors_Analysis.Rmd** - this notebook contains information about the data used and the code for all statistical analyses. 
- **2_Evolution-Multiple-Stressors_Figures.Rmd** - this notebook contains the code required to reproduce the figures (which were subsequently annotated in powerpoint). 
- **3_Evolution-Multiple-Stressors_SM.Rmd** - this notebook contains the code required to reproduce the supplementary figures and gifs used to make supplementary videos. 

There are also three folders: 

- **Data** - this folder contains the data used for all analyses in a csv file: *growth_data.csv* and will be populated by output data from the first notebook, which is required in the second and third notebooks. 
- **Scripts** - this folder contains a script with functions to calculate hedges effect sizes for interactions between stressors 
- **Images** - this folder contains an image (*growth_assay_table.png*) used to explain the data in the first notebook. 
- **gif** - this folder is populated by three gifs when the third notebook is run. These gifs were used to make the supplementary videos. 
