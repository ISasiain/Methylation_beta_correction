# METHYLATION DATA CORRECTION

#### Author: Iñaki Sasiain

### EXPLANATION ABOUT THE SCRIPT

### EXAMPLE 5000 CpGs

1. Splitting the data in training and validation and generating RData file with teh output

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/original_data; 
Rscript ../../scripts/get_data_to_analyse/extract_TV_from_Robject.r -d ./workspace_tcgaBrca_top5000.RData -b betaUnadjusted -p purityVector -s TRUE -v 20; 
```

2. Calculate the corrected values and regression parameters from the training data (it need 21 mins)

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/pop_regressions; 
Rscript ../../scripts/calculate_regs/new_purity_corrector.r -c 7 -b ../original_data/betas_training.RData -p ../original_data/purity_training.RData;
```

3. Use the regression parameters determined from the training data to estimate purities from validation data (it need 2 mins)

```bash
cd ~/Methylation/adjustBetas/01_5000_CpG/output;
Rscript ../../scripts/calculate_purity/run_all_validation.r -c 7 -d ../pop_regressions -b ../original_data/betas_validation.RData -o estimated_purity_5000CpG;
```
4. Analyse output data to produce plots. 

- An R script was used.

