# MAdaSub

Source files for the Metropolized Adaptive Subspace (MAdaSub) algorithm introduced in  

Staerk, C., Kateri, M., & Ntzoufras, I. (2021). A Metropolized adaptive subspace algorithm for high-dimensional Bayesian variable selection.

## Instructions

Please make sure that a version of R >= 3.6.0 is installed.

In order to replicate the results of the paper, please start R, 
set the working directory (setwd(..)) to the unzipped folder
and run the command: 

source("Master_file.R")

You have the options to replicate the analysis of the illustrative simulated example (Section 5.1), the simulation study (Section 5.2), the Tecator data application (Section 6.1), the PCR data application (Section 6.2) and the Leukemia data application (Section 6.3). 

The file "MAdaSub_LOAD.R" includes functions for running the serial MAdaSub algorithm ("MAdaSub") as well as its parallel version ("MAdaSub_parallel"). These functions can be used to examine further simulated or real data examples with MAdaSub. For instructions, please see the beginning of the file "MAdaSub_LOAD.R" (with description of input and output format).

For the PCR data application, please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)
from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip"
and extract the files Xgene.txt, gene_id.txt & Y3.txt from the zip file 77-5.Song.zip to the subfolder "Files".

Please note that the computation time for the full simulation study is quite long (considering 100 simulated datasets for each setting and considering the g-prior using the function "zlm" from the R-package "BMS", not optimized for computational speed). To obtain faster results, you also have the option to choose a lower number of simulated data examples. For the real data applications you also have the option to choose a smaller number of iterations as well as a smaller numbers of rounds and chains for the parallel MAdaSub algorithm. 

