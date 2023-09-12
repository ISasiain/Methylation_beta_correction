
# METHYLATION DATA CORRECTION

#### Author: Iñaki Sasiain

### SCRIPTS

### EXAMPLE 5000 CpGs

1. Splitting the data in training and validation and generating RData file with teh output

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/original_data; 
Rscript ../../scripts/get_data_to_analyse/preprocessing data.r -s TRUE-B ./workspace_tcgaBrca_top5000.RData -b betaUnadjusted -p purityVector -S TRUE -v 20; 
```

2. Calculate the corrected values and regression parameters from the training data (it need 21 mins)

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/pop_regressions; 
Rscript ../../scripts/calculate_regs/new_purity_corrector.r -c 7 -b ../original_data/betas_training.RData -p ../original_data/purity_training.RData;
```

3. Use the regression parameters determined from the training data to estimate purities from validation data (it needs 2 mins)

The smoothening and the correction of the correction has also been included

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/output;

#Calculating the results with the correction method
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regressions -b ../original_data/betas_validation.RData -o corr_estimated_purity_5000CpG -a 0.75 -r 0.6 -s 0.25 -p 5;

#Calulating the results without the correction method (The line that corrects the coverage was commented; line 45 of purity_coverage.r)
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regressions -b ../original_data/betas_validation.RData -o uncorr_estimated_purity_5000CpG -a 0.75 -r 0.6 -s 0.25 -p 5; NOT REPEATED

```
4. Analyse output data to produce plots. 

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/plots;

#The corrected data
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/corr_estimated_purity_5000CpG.RData -a ../original_data/purity_validation.RData -c ../output/corr_estimated_purity_5000CpG.used_cpgs.RData  -b ../original_data/betas_validation.RData  -o 5k_corr;

#The uncorrected data
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/uncorr_estimated_purity_5000CpG.RData -a ../original_data/purity_validation.RData -o 5k_uncorr;

#The corrected and smoothened data
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/corr_smooth_estimated_purity_5000CpG.RData -a ../original_data/purity_validation.RData -o 5k_smooth_corr;
```

5. Analyse the signal overestimation in the low purity region (this was done with the smootehned data);
 - Splitting the original validation dataset regarding their actual purity values
```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/original_data;
mkdir purity_splitted;
cd purity_splitted;
R;
```
-Split the dataset using R
```R
#Loading the data
betas <- readRDS("../betas_validation.RData");
purities <- readRDS("../purity_validation.RData");

#Creating a list with the 8 purity intrevals to split the data into (The names are 1-Purity, not purities)
list_to_split <- list(
    from_1_to_0.7 = c(0, 0.3),
    from_0.7_to_0.6 = c(0.3, 0.4),
    from_0.6_to_0.5 = c(0.4, 0.5),
    from_0.5_to_0.4 = c(0.5, 0.6),
    from_0.4_to_0.35 = c(0.6, 0.65),
    from_0.35_to_0.3 = c(0.65, 0.7),
    from_0.3_to_0.2 = c(0.7, 0.8),
    from_0.2_to_0 = c(0.8, 1)
)

#Iterate through the list to create a vector with the sample names included in each interval
for (interval in names(list_to_split)) {
    
    assign(paste("pur_",interval,sep=""),purities[purities >= list_to_split[[interval]][1] & purities < list_to_split[[interval]][2]])
    assign(paste("betas_",interval,sep=""), betas[,names(purities[purities >= list_to_split[[interval]][1] & purities < list_to_split[[interval]][2]])])

    saveRDS(get(paste("betas_",interval,sep="")), file=paste("betas_",interval,".RData",sep=""))
    saveRDS(get(paste("pur_",interval,sep="")), file=paste("pur_",interval,".RData",sep=""))
}
```

- Predicting purity per for each purity group. Creating a plot with all the coverage plots of the samples included into that actual purity interval. The purity_coverage function has to eb adapted to get the data to plot. Line 129 must be uncommented in order to get the required data.

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/plots;
mkdir analyse_overestimation/data_to_plot;
cd analyse_overestimation/data_to_plot;


#Getting the data tlo be plotted through a bash loop;
for int in $(ls ../../../original_data/purity_splitted/betas_*.RData); 
    do filename=$(echo ${int} | cut -d \/ -f 6 | sed 's/.RData//');
       echo ${filename};
       Rscript ../../../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../../../pop_regressions -b ${int} -o ${filename}; 
    done;

#Plotting the data using a bash loop and a R script
cd ../;

for my_file in $(ls ./data_to_plot/*);
    do plotname=$(echo ${my_file} | cut -d \/ -f 3 | sed 's/.RData//' | sed 's/betas_//' );
       Rscript plot_coverage.R -i ${my_file} -t ${plotname} -o ${plotname}
    done;

```

!!!! CHANGE THE LOCATION OF THE SCRIPTS !!!!!!


### 450.000 CpGs (All was run in Corsaire)

1. Transforming the purities from csv to Robject;

```bash
cd /home/Illumina/Iñaki_Sasiain/data;
conda activate meth_corr;
R;
```
Using R to transform the csv into an R object
```R
#Reading csv file
purities <- read.csv("TCGA_BRCA_tumorContent.csv");

#Creating a named vector from the csv file
purityVector <- purities$Cancer.DNA.fraction;
names(purityVector) <- purities$bcr_patient_barcode;

#Saving vector as a R object
save(purityVector, file="450k_CpGs_purities.RData")
```

2. Running a R script to remove CpGs from sexual chromosomes and split samples in training and validation datasets
```bash
cd /home/Illumina/Iñaki_Sasiain/data;
Rscript ../scripts/get_data_to_analyse/preprocessing data.r -s FALSE -B data450k_421368x630_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P 450k_CpGs_purities.RData -b betaOrig -p purityVector -S TRUE -v 20 -f TRUE -A object_450k_probesKeep.RData -a probesKeep -c chr
```

3. Calculating reference regressions from training dataset.
```bash
cd /home/Illumina/Iñaki_Sasiain/regressions;
nohup Rscript ../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ../data/betas_training.RData -p ../data/purity_training.RData &;
```

4. Estimate purities using the refernce regressions and the validation beta values.

```bash

#Using the same thresholds than in the first step
cd /home/Illumina/Iñaki_Sasiain/estimate_purities;
nohup Rscript ../scripts/calculate_purity/run_all_validation.r -c 35 -d ../regressions -b ../data/betas_validation.RData -o estimated_purity_450kCpG &;

#Using more strict thresholds (worse performance)
cd /home/Illumina/Iñaki_Sasiain/estimate_purity;
nohup Rscript ../scripts/calculate_purity/run_all_validation.r -c 25 -d ../regressions -b ../data/betas_validation.RData -o s1_estimated_purity_450kCpG -a 0.5 -r 0.4 -s 0.5 -p 5 &;

#Running it without filtering based on RSE and slope.
cd /home/Illumina/Iñaki_Sasiain/estimate_purity;
nohup Rscript ../scripts/calculate_purity/run_all_validation.r -c 30 -d ../regressions -b ../data/betas_validation.RData -o s2_estimated_purity_450kCpG -a 0.75 -r 0.6 -s 0.25 -p 5 &;
```


4. Plotting results;

```bash
#Using the same thresholds than in the first step
cd /home/Illumina/Iñaki_Sasiain/analyse_output;
Rscript ../scripts/analyse_output/analyse_output.r -e ../estimate_purity/estimated_purity_450kCpG.RData -a ../data/purity_validation.RData -p 450k_corr;

#Using more strict thresholds
cd /home/Illumina/Iñaki_Sasiain/analyse_output;
Rscript ../scripts/analyse_output/analyse_output.r -e ../estimate_purity/s1_estimated_purity_450kCpG.RData -a ../data/purity_validation.RData -p 450k_corr -o s1_result_analysis;

#Using different thesholds
cd /home/Illumina/Iñaki_Sasiain/analyse_output;
Rscript ../scripts/analyse_output/analyse_output.r -e ../estimate_purity/s2_estimated_purity_450kCpG.RData -a ../data/purity_validation.RData -p 450k_corr -o s2_result_analysis;
```

### Checking the effect of the beta estimation method (infinium I and II)

I HAVE CHANGED THE DIRECTORIES, ADAPT THE COMMANDS

1. Transforming the data into an R object

```R
methodology <- read.csv("infiniumInfo_421368CpGs.csv");
rownames(methodology) <- methodology$id;
saveRDS(methodology, file="infinium_info.RData")
```

2. Getting vectors of the CpGs quantified with each method in a list
```R
infinium_info <- readRDS("infinium_info.RData");

design_I <- rownames(infinium_info[infinium_info$designType=="I",]);
design_II <- rownames(infinium_info[infinium_info$designType=="II",]);

design_type <- list(
    "I" = design_I,
    "II" = design_II
)

saveRDS(design_type, file="designType_list.RData")
```

2. Creating training and validation datasets for each infinium method

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/original_data/infinium_splitted;
Rscript ../../../scripts/get_data_to_analyse/preprocessing_data.r -s TRUE -B ../workspace_tcgaBrca_top5000.RData -b betaUnadjusted -p purityVector -S TRUE -v 20 -d TRUE -D ../designType_list.RData;
```

3. Determining population regressions

```bash
cd ~/Methylation/adjustBetas/03_infinium/pop_regs;

#Creating refernce regression datasets for each infinium type
for inf in ../original_data/infinium_splitted/betas_training_*; 
    do Rscript ../../scripts/calculate_regs/new_purity_corrector.r -c 7 -b ${inf} -p ../original_data/infinium_splitted/purity_training.RData -o $(echo ${inf} | cut -d \/ -f 4 | sed 's/.RData//'); 
    done

#Creating different directories for the created output files and moving the files
mkdir infinium_I;
mkdir infinium_II;

mv betas_training_Ioutput_* ./infinium_I/;
mv betas_training_IIoutput_* ./infinium_II/;
```


4. Predicting purity
```bash
cd ~/Methylation/adjustBetas/03_infinium/output;

#Predicting from infinium I CpGs (correction + smoothening)
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regs/infinium_I -b ../original_data/infinium_splitted/betas_validation_I.RData -o corr_smooth_infI_5000CpG -a 0.75 -r 0.6 -s 0.25 -p 5;

#Predicting from infinium II CpGs (correction + smoothening)
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regs/infinium_II -b ../original_data/infinium_splitted/betas_validation_II.RData -o corr_smooth_infII_5000CpG -a 0.75 -r 0.6 -s 0.25 -p 5;
```

5. Creating plots per ecah infinium type
```bash
cd ~/Methylation/adjustBetas/03_infinium/plots;

#Producing plots for infinium I
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/corr_smooth_infI_5000CpG.RData -a ../original_data/infinium_splitted/purity_validation.RData -o 5k_corr_smooth_infI;

#Producing plots for infinium II
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/corr_smooth_infII_5000CpG.RData -a ../original_data/infinium_splitted/purity_validation.RData -o 5k_corr_smooth_infII;
```

### Estimating purity using different number of CpGs based on variability (Run in Corsaire)

1. Generating beta datasets with the 100, 250, 500, 1.000, 2.500, 5.000, 10.000, 50.000 and 100.000 most variable CpGs. The sexual chromosomes will be ignored

```bash
cd ~/Iñaki_Sasiain/04_CpG_nums/data;

#Defining cpg number list 
#cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);
cpg_list=(30000);

# Creating the training (80 %) and validation (20 %) datasets through a bash loop.
# The cpgs belonging to sexual chromosomes have been excluded
for num in ${cpg_list[@]}; 
    do mkdir ~/Iñaki_Sasiain/04_CpG_nums/data/cpgs_${num};
        cd ~/Iñaki_Sasiain/04_CpG_nums/data/cpgs_${num};
        Rscript ../../../scripts/get_data_to_analyse/preprocessing_data.r -s FALSE -B ../../../data/data450k_421368x630_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P ../../../data/450k_CpGs_purities.RData -b betaOrig -p purityVector -S TRUE -v 20 -A object_450k_probesKeep.RData -a probesKeep -c chr -N TRUE -n ${num};
    done;
```

2. Determining the regression for each beta dataset

```bash
cd ~/Iñaki_Sasiain/04_CpG_nums/calculate_regressions;

#Defining cpg number list 
cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);

#Calculating regressions for each cpg number using a bash loop
for num in ${cpg_list[@]}; 
    do mkdir /home/Illumina/Iñaki_Sasiain/04_CpG_nums/calculate_regressions/cpgs_${num};
        cd ~/Iñaki_Sasiain/04_CpG_nums/calculate_regressions/cpgs_${num};
        Rscript ../../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ~/Iñaki_Sasiain/04_CpG_nums/data/cpgs_${num}/betas_training.RData -p ~/Iñaki_Sasiain/04_CpG_nums/data/cpgs_${num}/purity_training.RData -o cpgs${num}; 
    done;

```

3. Estimating purity for each beta dataset

```bash
cd ~/Iñaki_Sasiain/04_CpG_nums/estimate_purity;

#Defining cpg number list 
cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);

#Estimating purity for each cpg number through a bash loop
for num in ${cpg_list[@]};
    do mkdir ~/Iñaki_Sasiain/04_CpG_nums/estimate_purity/cpgs_${num};
        cd ~/Iñaki_Sasiain/04_CpG_nums/estimate_purity/cpgs_${num};
        Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../calculate_regressions/cpgs_${num} -b ../../data/cpgs_${num}/betas_validation.RData -o cpgs_${num}.corr.smooth -a 0.75 -r 0.6 -s 0.25 -p 5;
    done;

```

4. Producing plots for each dataset
```bash
cd ~/Iñaki_Sasiain/04_CpG_nums/plots;

#Defining cpg number list 
cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);

# Producing plots for each cpg number through loop
for num in ${cpg_list[@]};
    do mkdir ~/Iñaki_Sasiain/04_CpG_nums/plots/cpgs_${num};
        cd ~/Iñaki_Sasiain/04_CpG_nums/plots/cpgs_${num};
        # IMPORTANT! Line 65 of analyse_output.r has been uncommented to adapt the sample name
        Rscript ../../../scripts/analyse_output/analyse_output.r -e ../../estimate_purity/cpgs_${num}/cpgs_${num}.corr.smooth.RData -a ../../data/cpgs_${num}/purity_validation.RData -c ../../estimate_purity/cpgs_${num}/cpgs_${num}.corr.smooth.used_cpgs.RData -b ../../data/cpgs_${num}/betas_validation.RData -o cpgs_${num};
    done;

# Comparing predictions
cd /home/Illumina/Iñaki_Sasiain/04_CpG_nums/plots;

paths_to_compare=$(for dir in ../estimate_purity/*; 
    do ls ${dir}/*.corr.smooth.RData | tr "\n" ",";
    done);

Rscript ../../scripts/analyse_output/compare_predictions.r -c ${paths_to_compare} -p ../data/cpgs_100/purity_validation.RData -o num_cpg_comparison;

# Plotting the percentages of the populations detcted in each cpg (using the slopes of the regressions)
cd /home/Illumina/Iñaki_Sasiain/04_CpG_nums/plots;

paths_to_compare=$(for dir in ../calculate_regressions/*; 
    do ls ${dir}/*output_reg.intercepts.RData | tr "\n" ",";
    done);

Rscript ../../scripts/analyse_output/evolution_pops_detected.r -f ${paths_to_compare[@]} -o num_cpgs;

```

### Estimating purity only with the CpGs used in all the sample's estimation (Corsaire 450k)

1. Copying original data and the regressions calculated for the 450k CpGs to the current directory
```bash
cd /home/Illumina/Iñaki_Sasiain/05_common_CpGs;

#Copying the directories
cp -r ../data .;
cp -r ../regressions .;
```

2. Estimating purity using only the common CpGs
```bash
cd /home/Illumina/Iñaki_Sasiain/05_common_CpGs/estimate_purity;

nohup Rscript  ../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../regressions -b ../data/betas_validation.RData -o s2_estimated_purity_450kCpG -a 0.75 -r 0.6 -s 0.25 -p 5 -i TRUE &;
```

3. Plotting the results
```bash
cd /home/Illumina/Iñaki_Sasiain/05_common_CpGs/plots;

Rscript ../../scripts/analyse_output/analyse_output.r -e ../estimate_purity/s2_estimated_purity_450kCpG.RData -a ../data/purity_validation.RData -c ../estimate_purity/s2_estimated_purity_450kCpG.used_cpgs.RData -b ../data/betas_validation.RData -o common_cpgs;

```
The results obtained were really unaccurate. The value obtained for the most of the samples correspons to the noise. May the common cpgs be responsible of the noise? In order to check that I will rerun the analysis only with the cpgs that are not common.

### Estimating purity only with the CpGs NOT used in all the sample's estimation (Corsaire 450k)

1. Copying original data and the regressions calculated for the 450k CpGs to the current directory

```bash
cd /home/Illumina/Iñaki_Sasiain/06_without_common_CpGs;

#Copying the directories
cp -r ../data .;
cp -r ../regressions .;
```

2. Estimate purity using only NOT COMMON CpGs. Line 300 was added and teh common_cpgs vector was changed for cpgs_to_use when recalculating the purity.

```bash
cd /home/Illumina/Iñaki_Sasiain/06_without_common_CpGs/estimate_purity;

nohup Rscript  ../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../regressions -b ../data/betas_validation.RData -o s2_estimated_purity_450kCpG -a 0.75 -r 0.6 -s 0.25 -p 5 -i TRUE &;
```

### Optimizing parameters (Run in corsaire)

1. Copying the data and regressions calculated to new directory (using 30.000 CpGs)

```bash
cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/data;
cp ../../04_CpG_nums/data/cpgs_30000/* .;

cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/regressions;
cp ../../04_CpG_nums/calculate_regressions/cpgs_30000/* .;
```

2. Estimating purity using different alpha values; 
```bash
cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/alpha_purities;

#Defining a list of alphas
alphas=(0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95)

for a in ${alphas[@]};
    do mkdir /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/alpha_purities/alpha_${a};
        cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/alpha_purities/alpha_${a};
        Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../regressions -b ../../data/betas_validation.RData -o alpha_${a}.corr.smooth -a ${a} -r 0.6 -s 0.25 -p 5;
    done;
```

3. Estimating purity using different slope threshold values; 
```bash
cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/slope_purities;

#Defining a list of slopes
slopes=(0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8)

for s in ${slopes[@]};
    do mkdir /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/slope_purities/slope_${s};
        cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/slope_purities/slope_${s};
        Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../regressions -b ../../data/betas_validation.RData -o slope_${s}.corr.smooth -a 0.75 -r 0.6 -s ${s} -p 5;
    done;
```

4. Estimating purity using different max RSE thresholds; 
```bash
cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/rse_purities;

#Defining a list of max RSE thresholds
rse=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 10000)

for r in ${rse[@]};
    do mkdir /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/slope_purities/rse_${r};
        cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/slope_purities/rse_${r};
        Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../regressions -b ../../data/betas_validation.RData -o rse_${r}.corr.smooth -a 0.75 -r ${r} -s 0.25 -p 5;
    done;
```

5. Comparing predictions
```bash
cd /home/Illumina/Iñaki_Sasiain/07_Estimating_parameters/comparing_results/;

# Comparing predictions made with different alphas
paths_to_compare=$(for dir in ../alpha_purities/*; 
    do ls ${dir}/*.corr.smooth.RData | tr "\n" ",";
    done);

Rscript ../../scripts/analyse_output/compare_predictions.r -c ${paths_to_compare} -p ../data/purity_validation.RData -o alpha_comparison;

# Comparing predictions made with different slope thresholds
paths_to_compare=$(for dir in ../slope_purities/*; 
    do ls ${dir}/*.corr.smooth.RData | tr "\n" ",";
    done);

Rscript ../../scripts/analyse_output/compare_predictions.r -c ${paths_to_compare} -p ../data/purity_validation.RData -o slope_comparison;

# Comparing predictions made with different RSE thresholds
paths_to_compare=$(for dir in ../rse_purities/*; 
    do ls ${dir}/*.corr.smooth.RData | tr "\n" ",";
    done);

Rscript ../../scripts/analyse_output/compare_predictions.r -c ${paths_to_compare} -p ../data/purity_validation.RData -o rse_comparison;
```

### Checking performance and CpG numbers to use through 6-fold-cross-validation (Run in corsaire)

1. Generating test and training datasets for each CpG number using cross validation;

```bash
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/data;

#Defining cpg number list 
cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);

# Producing plots for each cpg number through loop
for num in ${cpg_list[@]};
    do mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/data/cpgs_${num};
        cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/data/cpgs_${num};

        Rscript ../../../scripts/get_data_to_analyse/split_cross_validation.r -s FALSE -B ../../../data/data450k_421368x630_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P ../../../data/450k_CpGs_purities.RData -b betaOrig -u purityVector -S FALSE -C TRUE -k 6 -A object_450k_probesKeep.RData -a probesKeep -c chr -N TRUE -n ${num};
    done;
```


2. Generating the regressions for included cpg  number and fold

```bash
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/calculate_regressions;

nohup bash -c '
#Defining cpg number list 
cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);

#Creating the regerssion datasets using nohup. The process takes a very long time!
for num in ${cpg_list[@]}; 
    do mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/calculate_regressions/cpgs_${num};
       cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/calculate_regressions;
       for dir in $(ls ../data/cpgs_${num}/*_BetasTraining.RData); 
          do fold=$(echo ${dir} | cut -d \/ -f 4 | cut -d _ -f 1);
          mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/calculate_regressions/cpgs_${num}/${fold};
          cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/calculate_regressions/cpgs_${num}/${fold};
          Rscript ../../../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b /home/Illumina/Iñaki_Sasiain/08_Cross_validation/data/cpgs_${num}/${fold}_BetasTraining.RData -p /home/Illumina/Iñaki_Sasiain/08_Cross_validation/data/cpgs_${num}/${fold}_PurityTraining.RData -o ${fold};
    done;
done;
'
```

3. Estimating purity for cpg number and fold

```bash
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity;

#Running the scripts using nohup
nohup bash -c '

#Defining cpg number list 
cpg_list=(100 250 500 1000 2500 5000 10000 20000 30000 40000 50000 75000 100000 200000 421368);
#Determining purity for each cpg number and fold
for num in ${cpg_list[@]}; 
    do cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity;
       mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/cpgs_${num};
       for dir in $(ls ../data/cpgs_${num}/*_BetasTraining.RData); 
          do fold=$(echo ${dir} | cut -d \/ -f 4 | cut -d _ -f 1);
             mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/cpgs_${num}/${fold};
             cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/cpgs_${num}/${fold};
             Rscript ../../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../../calculate_regressions/cpgs_${num}/${fold} -b ../../../data/cpgs_${num}/${fold}_BetasTest.RData  -o PredPurity_${fold}_cpgs${num} -a 0.75 -s 0.25 -p 5;
    done;
done;
'
```

4. Estimating purity for alpha value and fold

# The regressions that were estimated to calculate the cross validation for the cpg number will be used.


# Eatimating purities. Using 30.000 as the cpgs number for the predictions
```bash
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity;

#Running the scripts using nohup
nohup bash -c '

#Defining alpha list 
alphas=(0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95);
#Determining purity for each alpha and fold
for num in ${alphas[@]}; 
    do cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity;
       mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/alpha_${num};
       for dir in $(ls ../data/cpgs_30000/*_BetasTraining.RData); #Using 30000 cpgs in the prediction
          do fold=$(echo ${dir} | cut -d \/ -f 4 | cut -d _ -f 1);
             mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/alpha_${num}/${fold};
             cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/alpha_${num}/${fold};
             Rscript ../../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../../calculate_regressions/cpgs_30000/${fold} -b ../../../data/cpgs_30000/${fold}_BetasTest.RData -o PredPurity_${fold}_alpha${num} -a ${num} -s 0.25 -p 5;
    done;
done;
'
```

5. Estimating purity for slope threshold and fold

# The regressions that were estimated to calculate the cross validation for the cpg number will be used.


# Eatimating purities. Using 30.000 as the cpgs number for the predictions
```bash
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity;

#Running the scripts using nohup
nohup bash -c '

#Defining slope threshold list 
slopes=(0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8)
#Determining purity for each slope threshold and fold
for num in ${slopes[@]}; 
    do cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity;
       mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/slope_${num};
       for dir in $(ls ../data/cpgs_30000/*_BetasTraining.RData); #Using 30000 cpgs in the prediction
          do fold=$(echo ${dir} | cut -d \/ -f 4 | cut -d _ -f 1);
             mkdir /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/slope_${num}/${fold};
             cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity/slope_${num}/${fold};
             Rscript ../../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../../calculate_regressions/cpgs_30000/${fold} -b ../../../data/cpgs_30000/${fold}_BetasTest.RData -o PredPurity_${fold}_slope${num} -a 0.75 -s ${num} -p 5;
    done;
done;
'
```

4. Comparing results. 

```bash

# Creating a R object with the actual purity 
cd /home/Illumina/Iñaki_Sasiain/data;
Rscript ../scripts/get_data_to_analyse/preprocessing:data.r -s FALSE -B data450k_421368x630_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P 450k_CpGs_purities.RData -b betaOrig -p purityVector -S FALSE -f FALSE -N FALSE; 


 #                             #
 # COMPARE RESULTS: CPG NUMBER #
 #                             #

# Define variables required to run the script
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/analyse_output;
list_of_folds=$(find ../estimate_purity/cpgs_100 -mindepth 1 -maxdepth 1 -type d | cut -d \/ -f 4 | tr "\n" "," | sed 's/,$//');
list_of_cpgs=$(find ../estimate_purity/cpgs* -maxdepth 0 | cut -d \/ -f 3 | tr "\n" "," | sed 's/,$//');

path="/home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity";

#Potting the results
Rscript /home/Illumina/Iñaki_Sasiain/scripts/analyse_output/compare_cross_validation.r -f ${list_of_folds} -c ${list_of_preds} -d $path -p /home/Illumina/Iñaki_Sasiain/data/purity.RData -o cpg_num_CrossVal; 


 #                        #
 # COMPARE RESULTS: ALPHA #
 #                        #

# Define variables required to run the script
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/analyse_output;
list_of_folds=$(find ../estimate_purity/cpgs_100 -mindepth 1 -maxdepth 1 -type d | cut -d \/ -f 4 | tr "\n" "," | sed 's/,$//');
list_of_alphas=$(find ../estimate_purity/alpha* -maxdepth 0 | cut -d \/ -f 3 | tr "\n" "," | sed 's/,$//');

path="/home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity";

#Potting the results
Rscript /home/Illumina/Iñaki_Sasiain/scripts/analyse_output/compare_cross_validation.r -f ${list_of_folds} -c ${list_of_alphas} -d ${path} -p /home/Illumina/Iñaki_Sasiain/data/purity.RData -o alpha_CrossVal;


 #                        #
 # COMPARE RESULTS: SLOPE #
 #                        #

# Define variables required to run the script
cd /home/Illumina/Iñaki_Sasiain/08_Cross_validation/analyse_output;
list_of_folds=$(find ../estimate_purity/cpgs_100 -mindepth 1 -maxdepth 1 -type d | cut -d \/ -f 4 | tr "\n" "," | sed 's/,$//');
list_of_slopes=$(find ../estimate_purity/slope* -maxdepth 0 | cut -d \/ -f 3 | tr "\n" "," | sed 's/,$//');

path="/home/Illumina/Iñaki_Sasiain/08_Cross_validation/estimate_purity";

#Potting the results
Rscript /home/Illumina/Iñaki_Sasiain/scripts/analyse_output/compare_cross_validation.r -f ${list_of_folds} -c ${list_of_slopes} -d ${path} -p /home/Illumina/Iñaki_Sasiain/data/purity.RData -o slope_CrossVal; 
```


### Testing the tool with new data (Run in Corsaire)

#### Using GSE148748 as new test data

1. Getting and preprocessing GSE148748 data

```bash
# Creating new directory to store the new data
mkdir /home/Illumina/Iñaki_Sasiain/data/GSE148748_data;
cd /home/Illumina/Iñaki_Sasiain/data/GSE148748_data;

# Getting and unzipping the data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148748/matrix/GSE148748_series_matrix.txt.gz;
gunzip GSE148748/matrix/GSE148748_series_matrix.txt.gz;

# Preprocessing the downloaded file in R
```R
# Converting the tsv in a R object and removing headers.
GSE148748 <- read.csv("GSE148748_series_matrix.txt", sep="\t", skip=56, na.strings=c("", NA), head=1);

# Setting the sample and cpg ids as colnames and rownames
rownames(GSE148748) <- GSE148748$"ID_REF";
GSE148748 <- GSE148748[,-which(names(GSE148748)=="ID_REF")];

# Removing cpgs without beta values
GSE148748 <- GSE148748[-which(rowSums(is.na(GSE148748))==ncol(GSE148748)),];

#Saving the dataframe as an R object
saveRDS(GSE148748, file="GSE148748_betas.RData");
```

2. Getting test and training data to run the analysis

```bash
#Using all the cpgs

# Getting training dataset (TNBC TCGA data)
cd /home/Illumina/Iñaki_Sasiain/09_TNBC_final/data/training;
cp ../../../data/betas.RData .; # Complete betas dataset (450K CpG)
cp ../../../data/purity.RData .; # Complete purity dataset

# Getting test dataset (GSE148748)
cd /home/Illumina/Iñaki_Sasiain/09_TNBC_final/data/test;
cp ../../../data/GSE148748_data/GSE148748_betas.RData .; # Complete betas dataset

# Getting test dataset. Using only the 30.000 most variant CpGs of the refernce dataset
cd /home/Illumina/Iñaki_Sasiain/09_TNBC_final/data/test_30000;

Rscript ../../../scripts/get_data_to_analyse/get_most_variables_cpgs.r -r ../training/betas.RData -a ../test/GSE148748_betas.RData -n 30000 -l TRUE -p TNBC;
```

3. Using the whole TNBC samples of TCGA-Breast cancer to determint the reference regressions

```bash
cd /home/Illumina/Iñaki_Sasiain/09_TNBC_final/regressions;
Rscript ../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ../data/training/betas.RData -p ../data/training/purity.RData -o ref_reg_TNBC;
```

4. Estimating purity for GSE148748 based on the reference regressions

```bash
# Using all the CpGs to estimate purity
cd /home/Illumina/Iñaki_Sasiain/09_TNBC_final/estimating_purity/all_cpgs;
Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../regressions/ -b ../../data/test/GSE148748_betas.RData -o GSE148748_est_pur -a 0.75 -s 0.25 -p 5;

#Using only the 30000 most variable reference of the reference dataset
cd /home/Illumina/Iñaki_Sasiain/09_TNBC_final/estimating_purity/30000_cpgs;
Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 35 -d ../../regressions/ -b ../../data/test_30000/TNBC_most_variable_CpGs.RData -o GSE148748_est_pur -a 0.75 -s 0.25 -p 5;
```

#### Using LUAD data from TGCA for training and test

1. Getting the data. Obtained from the group and reformatted.
 
```R
#The LUAD_LUSC_purity.RData file was first reformatted to meet the required format
load("LUAD_LUSC_purity.RData")

purity_LUAD <- setNames(object=as.vector(purity_LUAD[,3]), as.vector(purity_LUAD[,1]))
purity_LUSC <- setNames(object=as.vector(purity_LUSC[,3]), as.vector(purity_LUSC[,1]))

save(purity_LUAD, purity_LUSC, file="LUAD_LUSC_purity.RData")
```

``` bash
# Getting complete dataset (without splitting in training and test)
cd /home/Illumina/Iñaki_Sasiain/10_LUAC_final/data/full_data;
Rscript ../../../scripts/get_data_to_analyse/preprocessing_data.r -s FALSE -B /home/Illumina/Iñaki_Sasiain/data/LUAD_data/LUAD_data450k_421368x418_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P /home/Illumina/Iñaki_Sasiain/data/LUAD_data/LUAD_LUSC_purity.RData -b betaOrig -p purity_LUAD -S FALSE -f FALSE -N FALSE; 

# Getting splitted dataset. 80% training 20% test.
cd /home/Illumina/Iñaki_Sasiain/10_LUAC_final/data/training_test;
Rscript ../../../scripts/get_data_to_analyse/preprocessing_data.r -s FALSE -S TRUE -v 20 -B /home/Illumina/Iñaki_Sasiain/data/LUSC_data/LUSC_data450k_421368x333_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P /home/Illumina/Iñaki_Sasiain/data/LUSC_data/LUAD_LUSC_purity.RData -b betaOrig -p purity_LUSC -f FALSE -N FALSE; 
```

2. Calculating regressions. TO RUN!!!!!!!

```bash
# Getting the regerssions for the complete dataset
cd /home/Illumina/Iñaki_Sasiain/10_LUAC_final/regressions/full_data;
Rscript ../../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ../../data/full_data/betas.RData -p ../../data/full_data/purity.RData -o ref_reg_LUAD;

#Getting the regressions for only with the splitted validation dataset
cd /home/Illumina/Iñaki_Sasiain/10_LUAC_final/regressions/training_test;
Rscript ../../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ../data/training_test/betas.RData -p ../data/training_test/purity.RData -o splitted_reg_LUAD;
```


#### Using LUSC data from TCGA for training and test

1. Getting the data. Obtained from the group and reformatted.

``` bash
# Getting complete dataset (without splitting in training and test)
cd /home/Illumina/Iñaki_Sasiain/11_LUSC_final/data/full_data;
Rscript ../../../scripts/get_data_to_analyse/preprocessing_data.r -s FALSE -B /home/Illumina/Iñaki_Sasiain/data/LUSC_data/LUSC_data450k_421368x333_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P /home/Illumina/Iñaki_Sasiain/data/LUSC_data/LUAD_LUSC_purity.RData -b betaOrig -p purity_LUAD -S FALSE -f FALSE -N FALSE; 

# Getting splitted dataset. 80% training 20% test.
cd /home/Illumina/Iñaki_Sasiain/11_LUSC_final/data/training_test;
Rscript ../../../scripts/get_data_to_analyse/preprocessing_data.r -s FALSE -S TRUE -v 20 -B /home/Illumina/Iñaki_Sasiain/data/LUSC_data/LUSC_data450k_421368x333_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData -P /home/Illumina/Iñaki_Sasiain/data/LUSC_data/LUAD_LUSC_purity.RData -b betaOrig -p purity_LUSC -f FALSE -N FALSE;
```

2. Calculating regressions. TO RUN!!!!!

```bash
# Getting the regerssions for the complete dataset
cd /home/Illumina/Iñaki_Sasiain/11_LUSC_final/regressions/full_data;
Rscript ../../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ../data/full_data/betas.RData -p ../data/full_data/purity.RData -o ref_reg_LUSC;

#Getting the regressions for only with the splitted validation dataset
cd /home/Illumina/Iñaki_Sasiain/11_LUSC_final/regressions/training_test;
Rscript ../../../scripts/calculate_regs/new_purity_corrector.r -c 35 -b ../data/training_test/betas.RData -p ../data/training_test/purity.RData -o splitted_reg_LUSC;
```