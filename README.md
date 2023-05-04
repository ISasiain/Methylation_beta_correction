# METHYLATION DATA CORRECTION

#### Author: Iñaki Sasiain

### SCRIPTS

### EXAMPLE 5000 CpGs

1. Splitting the data in training and validation and generating RData file with teh output

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/original_data; 
Rscript ../../scripts/get_data_to_analyse/preprocessing data.r -d ./workspace_tcgaBrca_top5000.RData -b betaUnadjusted -p purityVector -s TRUE -v 20; 
```

2. Calculate the corrected values and regression parameters from the training data (it need 21 mins)

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/pop_regressions; 
Rscript ../../scripts/calculate_regs/new_purity_corrector.r -c 7 -b ../original_data/betas_training.RData -p ../original_data/purity_training.RData;
```

3. Use the regression parameters determined from the training data to estimate purities from validation data (it needs 2 mins)

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/output;

#Calculating the results with the correction method
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regressions -b ../original_data/betas_validation.RData -o corr_estimated_purity_5000CpG;

#Calulating the results without the correction method (The line that corrects the coverage was commented; line 45 of purity_coverage.r)
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regressions -b ../original_data/betas_validation.RData -o uncorr_estimated_purity_5000CpG;
```
4. Analyse output data to produce plots. 

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/plots;

#The corrected data
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/corr_estimated_purity_5000CpG.RData -a ../original_data/purity_validation.RData -o 5k_corr;

#The uncorrected data
Rscript ../../scripts/analyse_output/analyse_output.r -e ../output/uncorr_estimated_purity_5000CpG.RData -a ../original_data/purity_validation.RData -o 5k_uncorr;
```

5. Analyse the signal overestimation in the low purity region;
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

#Creating a list with the 8 purity intrevals to split the data into
list_to_split <- list(
    one_of_eight=c(0, 0.2),
    two_of_eight=c(0.2, 0.4),
    three_of_eight=c(0.4, 0.6),
    four_of_eight=c(0.6, 0.8),
    five_of_eight=c(0.8, 1)
)

#Iterate through the list to create a vector with the sample names included in each interval
for (interval in names(list_to_split)) {
    
    assign(paste("pur_",interval,sep=""),purities[purities >= list_to_split[[interval]][1] & purities < list_to_split[[interval]][2]])
    assign(paste("betas_",interval,sep=""), betas[,names(purities[purities >= list_to_split[[interval]][1] & purities < list_to_split[[interval]][2]])])

    saveRDS(get(paste("betas_",interval,sep="")), file=paste("betas_",interval,".RData",sep=""))
    saveRDS(get(paste("pur_",interval,sep="")), file=paste("pur_",interval,".RData",sep=""))
}
```

- Predicting purity per for each purity group. Creating a plot with all the coverage plots of the samples included into that actual purity interval.

```bash
Rscript ../../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../../pop_regressions -b betas_one_of_eight.RData -o one_of_eight;
```


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
cd /home/Illumina/Iñaki_Sasiain/estimate_purities;
nohup Rscript ../scripts/calculate_purity/run_all_validation.r -c 35 -d ../regressions -b ../data/betas_validation.RData -o estimated_purity_450kCpG &;
```

4. Plotting results;

```bash
cd /home/Illumina/Iñaki_Sasiain/;
Rscript ../scripts/analyse_output/analyse_output.r -e ../estimate_purity/estimated_purity_450kCpG.RData -a ../data/purity_validation.RData -p 450k_corr;
```