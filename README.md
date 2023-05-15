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

### Estimating purity using different number of CpGs based on variability

1. Generating beta datasets with the 100, 250, 500, 1.000, 2.500, 5.000, 10.000, 50.000 and 100.000 most variable CpGs

```bash
cd ~/Methylation/adjustBetas/04_CpG_nums/original_data;

Rscript ../../scripts/get_data_to_analyse/preprocessing_data.r -s TRUE -B ../../01_5000_CpG/original_data/workspace_tcgaBrca_top5000.RData -b betaUnadjusted -p purityVector -S TRUE -v 20 -N TRUE -n 100
```