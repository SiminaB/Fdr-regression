#README for FDR-regression repository

This folder contains the code to fully reproduce the simulations and data analysis in the paper
"A direct approach to estimating false discovery rates conditional on covariates" by SM Boca and JT Leek.

Before running the code, please make sure you install the packages swfdr and FDRreg, available at:
https://github.com/leekgroup/swfdr    
https://github.com/jgscott/FDRreg

The .pdf files can be generated from the .Rmd/.Rnw files by opening them in RStudio version 0.99.903 or later and pressing "Compile PDF". In order for this to properly work, first go to Tools -> Global Options -> Sweave and select "knitr" (instead of "Sweave") for "Weave Rnw files using:". For "Typeset LaTeX into PDF using:" the "pdfLaTeX" option was selected. A LaTeX distribution such as MikTeX is required.

Some of the code is parallelized in order to substantially speed it up. Please note that this is only checked on a Windows machine and adjustments may need to be made for other operating systems and for machines with different numbers of cores.

Note that the simulation code and the bootstrapping code can be particularly computationally intensitive and some of the files may take several hours to run. 
Some of the output files are not included, due to reposity size limitations.

The directory structure is as follows:

###functions.R
The file **functions.R** contains a number of R functions used in some of the Rnw and Rmd files below.

##BMI GIANT GWAS meta-analysis
This directory contains the code to perform the data analysis and generate Figures 1, 2, and S5 which are saved in the "Figures" subdirectory.
The data must first be downloaded from https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_BMI
It is not included here, due to repository size limitations.

The file **1.read_in_GIANT.Rmd** reads in the data in tsv format and creates the **BMI_GIANT_GWAS.RData file**. Running it also generates the file **1.read_in_GIANT.pdf**.

The file **2.make_Figure_1.Rmd** loads BMI_GIANT_GWAS.RData and generates Figures 1 and S5 from the paper, saved in Figures/Fig3-1.pdf, Figures/FigS5-1.pdf, Figures/FigS5-2.pdf, as well as the numbers used in Table 2. Running it also generates the file 2.make_Figure_1.pdf.

The file **3.run_analysis.Rmd** loads BMI_GIANT_GWAS.RData and runs our approach, saving the results to **BMI_GIANT_GWAS_results_logistic.RData**. Running it also generates the file **3.run_analysis.pdf**.

The file **4.run_analysis_Scott.Rmd** loads BMI_GIANT_GWAS.RData and runs the Scott approach with the empirical null, saving the results to **BMI_GIANT_GWAS_results_Scott.RData**. Running it also generates the file **4.run_analysis_Scott.pdf**.

The file **4.run_analysis_Scott_theoretical.Rmd** loads BMI_GIANT_GWAS.RData and runs the Scott approach with the theoretical null, saving the results to **BMI_GIANT_GWAS_results_Scott_theoretical.RData**. Running it also generates the file **4.run_analysis_Scott_theoretical.pdf**.

The file **4b.bootstrap_CIs.Rmd** loads BMI_GIANT_GWAS.RData, BMI_GIANT_GWAS_results_logistic.RData and BMI_GIANT_GWAS_results_Scott.RData and runs 100 bootstrap iterations to generate a 90% confidence interval for the true proportion of nulls with our approach. The bootstrap results are saved to **BMI_GIANT_GWAS_bootstrap_all_logistic.RData** (not included here) and the bootstrap CI upper and lower limits to **BMI_GIANT_GWAS_bootstrap_CIs_logistic.RData**. Running it also generates **4b.bootstrap_CIs.pdf**. Note that running this can take 12-24 hours on a PC.

The file **5.make_Figure_2.Rmd** loads BMI_GIANT_GWAS.RData, BMI_GIANT_GWAS_results_logistic.RData, BMI_GIANT_GWAS_results_Scott.RData,  BMI_GIANT_GWAS_results_Scott_theoretical.RData, and **BMI_GIANT_GWAS_bootstrap_CIs_logistic.RData** generates Figure 2 from the paper, saved in **Figures/Fig2-1.pdf**. Running it also generates the file **5.make_Figure_2.pdf**.

##Simulations - independent
This directory contains the code to run the simulations for the independent test statistics, m=1,000 feature case, and generate Table 3 and Figures S1-S2. 
Before running it, should first create the subdirectories that start with alt*, which represent the different mechanisms for generating p-values under the alternative distribution.
Some of the output files are not included here due to repository size limitations.

The file **1.run_simulations.Rnw** simulates data and creates the files **alt\*/simResults_\*.RData**. Running it also generates the file **1.run_simulations.pdf**.

The file **2.estimate_pi0x_thresh.Rnw** inputs the alt\*/simResults\*.RData files and generates the **alt\*/simResults_pi0x_thresh\*.RData** files. Running it also generates the file **2.estimate_pi0x_thresh.pdf**.

The file **3.estimate_pi0x_Scott.Rnw** inputs the alt\*/simResults\*.RData files and generates the **simResults_pi0x_Scott\*.RData** files. Running it also generates the file **3.estimate_pi0x_Scott.pdf**.

The file **3.estimate_pi0x_Scott_emp.Rnw** inputs the alt\*/simResults\*.RData files and generates the **simResults_pi0x_Scott_emp\*.RData** files. Running it also generates the file **3.estimate_pi0x_Scott_emp.pdf**.

The file **4.make_FDR_TPR_tables.Rnw** inputs all the files in alt\* and generates the **alt\*/FDR_TPR_sims.RData** files. Running it also generates the file **4.make_FDR_TPR_tables.pdf**.

The file **5.make_combined_FDR_TDR_table.Rnw** inputs alt\*/FDR_TPR_sims.RData** files and generates tables with all the simulation results. Running it also generates the file **5.make_combined_FDR_TPR_table.pdf**.

The file **5.make_combined_FDR_TDR_table_only_WS.Rnw** inputs alt\*/FDR_TPR_sims.RData** files and generates Table 3. Running it also generates the file **5.make_combined_FDR_TPR_table_only_WS.pdf**.

The file **6.make_Figure_comparing_means.Rnw** inputs alt_z_large/simResults* and alt_t_large/simResults* and generates Figures S1 and S2. Running it also generates the file **6.make_Figure_comparing_means.pdf**.

##Simulations - independent - 10000
This directory contains the code to run the simulations for the independent test statistics, m=10,000 feature case, and generate Table 4 and Figures S3-S4. 
Before running it, should first create the subdirectories that start with alt*, which represent the different mechanisms for generating p-values under the alternative distribution.

Some of the output files are not included here due to repository size limitations.

The file **1.run_simulations.Rnw** simulates data and creates the files **alt\*/simResults_\*.RData**. Running it also generates the file **1.run_simulations.pdf**.

The file **2.estimate_pi0x_thresh.Rnw** inputs the alt\*/simResults\*.RData files and generates the **alt\*/simResults_pi0x_thresh\*.RData** files. Running it also generates the file **2.estimate_pi0x_thresh.pdf**.

The file **3.estimate_pi0x_Scott.Rnw** inputs the alt\*/simResults\*.RData files and generates the **simResults_pi0x_Scott\*.RData** files. Running it also generates the file **3.estimate_pi0x_Scott.pdf**.

The file **3.estimate_pi0x_Scott_emp.Rnw** inputs the alt\*/simResults\*.RData files and generates the **simResults_pi0x_Scott_emp\*.RData** files. Running it also generates the file **3.estimate_pi0x_Scott_emp.pdf**.

The file **4.make_FDR_TPR_tables.Rnw** inputs all the files in alt\* and generates the **alt\*/FDR_TPR_sims.RData** files. Running it also generates the file **4.make_FDR_TPR_tables.pdf**.

The file **5.make_combined_FDR_TDR_table.Rnw** inputs alt\*/FDR_TPR_sims.RData** files and generates tables with all the simulation results. Running it also generates the file **5.make_combined_FDR_TPR_table.pdf**.

The file **5.make_combined_FDR_TDR_table_only_WS.Rnw** inputs alt\*/FDR_TPR_sims.RData** files and generates Table 4. Running it also generates the file **5.make_combined_FDR_TPR_table_only_WS.pdf**.

The file **6.make_Figure_comparing_means.Rnw** inputs alt_z_large/simResults* and alt_t_large/simResults* and generates Figures S3 and S4. Running it also generates the file **6.make_Figure_comparing_means.pdf**.

