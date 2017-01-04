#README for FDR-regression repository

This folder contains the code to fully reproduce the simulations and data analysis in the paper
"A regression framework for the proportion of true null hypotheses" by SM Boca and JT Leek.

Before running the code, please make sure you install the packages swfdr and FDRreg, available at:
https://github.com/leekgroup/swfdr    
https://github.com/jgscott/FDRreg

Some of the code is parallelized in order to substantially speed it up. Please note that this is only checked on a Windows machine and adjustments may need to be made for other operating systems and for machines with different numbers of cores.

Note that the simulation code can be particularly computationally intensitive and some of the files may take several hours to run. (Can compare date last modified for .Rnw files and .pdf files to obtain an estimate.)

The directory structure is as follows:

The file functions.R contains a number of R functions used in some of the Rnw files below.

##BMI GIANT GWAS meta-analysis
This directory contains the code to perform the data analysis and generate Figures 3 and 4, which are saved in the "Figures" subdirectory.
The data must first be downloaded from https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_BMI
It is not included here, due to repository size limitations.

The file 1.read_in_GIANT.Rmd reads in the data in tsv format and creates the BMI_GIANT_GWAS.RData file. Running it also generates the file 1.read_in_GIANT.pdf.

The file 2.make_Figure_3.Rmd loads BMI_GIANT_GWAS.RData and generates Figure 3 from the paper, saved in Figures/Fig3-1.pdf. Running it also generates the file 2.make_Figure_3.pdf.

The file 3.run_analysis.Rmd loads BMI_GIANT_GWAS.RData and runs our approach for estimating the proportion of true null hypotheses in a regression framework, saving the results to BMI_GIANT_GWAS_results.RData. Running it also generates the file 3.run_analysis.pdf.

The file 4.run_analysis_Scott.Rmd loads BMI_GIANT_GWAS.RData and runs the Scott approach for estimating the proportion of true null hypotheses in a regression framework, saving the results to BMI_GIANT_GWAS_results_Scott.RData. Running it also generates the file 3.run_analysis_Scott.pdf.

The file 5.make_Figure_4.Rmd loads BMI_GIANT_GWAS.RData, BMI_GIANT_GWAS_results.RData, and BMI_GIANT_GWAS_results_Scott.RData and generates Figure 4 from the paper, saved in Figures/Fig4-1.pdf. Running it also generates the file 5.make_Figure_4.pdf.

##Simulations - independent
This directory contains the code to run the simulations and generate Figures 1, 2, and S1, and Table S1. 

The file 1.run_simulations.Rnw simulates data and creates the files simResults_1.RData, simResults_2.RData, simResults_3.RData. Each of these files is around 70 MB. They are not included here due to repository size limitations. Running it also generates the file 1.run_simulations.pdf.

The file 2.estimate_pi0x_noThresh.Rnw inputs the simResults*.RData files and generates the simResults_pi0x_noThresh*.RData files. Running it also generates the file 2.estimate_pi0x_noThresh.pdf.

The file 2.estimate_pi0x_thresh.Rnw inputs the simResults*.RData files and generates the simResults_pi0x_thresh*.RData files. Running it also generates the file 2.estimate_pi0x_thresh.pdf.

The file 3.estimate_pi0x_Scott.Rnw inputs the simResults*.RData files and generates the simResults_pi0x_Scott*.RData files. Running it also generates the file 3.estimate_pi0x_Scott.pdf.

The file 4.make_Figure_1_S1.Rnw inputs the simResults_pi0x_noThresh*.RData and simResults_pi0x_Scott*.RData files and generates Figures 1 and S1 from the paper, saved in Figures/Fig1*-1.pdf and Figures/FigS1*-1.pdf. Running it also generates the file 4.make_Figure_1_S1.pdf.

The file 4.make_Figure_2.Rnw inputs the simResults_pi0x_thresh*.RData and simResults_pi0x_Scott*.RData files and generates Figure 2 from the paper, saved in Figures/Fig2*-1.pdf. Running it also generates the file 4.make_Figure_2.pdf.

The file 5.make_Table_S1.Rnw generates the results that are used in Table S1 in the paper. The table can be obtained from the generated file 4.make_Table_S1.pdf.

##Simulations - dependent
This directory contains the code to run the simulations and generate Figures S2 and S3. 

The file 1.run_simulations.Rnw simulates data and creates the files simResults_1.RData, simResults_2.RData, simResults_3.RData. Each of these files is around 70 MB. They are not included here due to repository size limitations. Running it also generates the file 1.run_simulations.pdf.

The file 2.estimate_pi0x_noThresh.Rnw inputs the simResults*.RData files and generates the simResults_pi0x_noThresh*.RData files. Running it also generates the file 2.estimate_pi0x_noThresh.pdf.

The file 2.estimate_pi0x_thresh.Rnw inputs the simResults*.RData files and generates the simResults_pi0x_thresh*.RData files. Running it also generates the file 2.estimate_pi0x_thresh.pdf.

The file 3.estimate_pi0x_Scott.Rnw inputs the simResults*.RData files and generates the simResults_pi0x_Scott*.RData files. Running it also generates the file 3.estimate_pi0x_Scott.pdf.

The file 4.make_Figure_S2.Rnw inputs the simResults_pi0x_noThresh*.RData and simResults_pi0x_Scott*.RData files and generates Figure S2 from the paper, saved in Figures/FigS2*-1.pdf. Running it also generates the file 4.make_Figure_S2.pdf.

The file 4.make_Figure_S3.Rnw inputs the simResults_pi0x_thresh*.RData and simResults_pi0x_Scott*.RData files and generates Figure S3 from the paper, saved in Figures/FigS3*-1.pdf. Running it also generates the file 4.make_Figure_S3.pdf.

