This file shows the usage of the complete beta correction pipeline for samples without known purity values. We recommend to run this scripts in a conda environment with the required packages, which are detailed in the spec.txt file available at the project's GItHUb page.

1. Generation of reference regressions
```bash

cd ./generating_regressions; 

#Running the reference regression calculator script using 0.05 as the variance threshold usinbg 6 cores
cd ../var_0.05;
Rscript ../../../scripts/calculate_regs/ref_regression_calculator.r -c 6 -b ../../data/reference_data/betas_ref.rds -p ../../data/reference_data/purity_ref.rds -o ./ -n BRCA_var0.05 -v 0.05;

#Running the reference regression calculator script using 0.05 as the variance threshold using 6 cores
cd ../var_0;
Rscript ../../../scripts/calculate_regs/ref_regression_calculator.r -c 6 -b ../../data/reference_data/betas_ref.rds -p ../../data/reference_data/purity_ref.rds -o ./ -n BRCA_var0 -v 0;
```

2. Sample purity estimation

```bash
cd ../../estimating_purities;

#Running the purity estimation script using 6 cores
Rscript ../../scripts/calculate_purity/purity_estimator.r -c 6 -a 0.7 -s 0.25 -p 4 -d ../generating_regressions/var_0.05/ -b ../data/data_to_correct/betas_toCorrect.rds -o example_estimated_purity;
```

3. Beta correction

```bash
#Correcting betas refitting the reference regressions
cd ../correcting_betas/refitting;

Rscript ../../../scripts/final_beta_correction/final_beta_correction.r -c 6 -F FALSE -o ./ -n example_refitting -r TRUE -P ../../data/reference_data/purity_ref.rds -B ../../data/reference_data/betas_ref.rds -b ../../data/data_to_correct/betas_toCorrect.rds -p ../../estimating_purities/example_estimated_purity.tsv

#Correcting betas without refitting the reference regressions
cd ../not_refitting;

Rscript ../../../scripts/final_beta_correction/final_beta_correction.r -c 6 -F FALSE -o ./ -n example_refitting -r FALSE -R ../../generating_regressions/var_0/ -b ../../data/data_to_correct/betas_toCorrect.rds -p ../../estimating_purities/example_estimated_purity.tsv;
``````
