# Fogli Iseppe et al. 2021 male and female human epicardial ventricular cardiomyocyte models and ML code.

The parent model developed by Yang and Clancy (Yang, P.-C. & Clancy, C. E. In silico Prediction of Sex-Based Differences in Human Susceptibility to Cardiac Ventricular Tachyarrhythmias. Front. Physiol. 3, 1–12 (2012)) was modified by:
- increasing the maximal transport rate of the Na+-Ca2+ exchanger (NCX) by 15% in the female model
- removing the original female-to-male difference in SERCA formulation
- removing the original female-to-male difference in Na+/K+-ATPase formulation
The updated model better recapitulates observed functional sex differences in Ca2+ handling.

_______________________________________________________________________________________________________________

## Contents:

- *readme.txt*		 		this file

### "AFI_model" folder

- *run_AFI_model.py*			runs the simulations with the baseline male and female models
- *AFI_ODEfile.cpp*			source file containing ODEs for the Fogli Iseppe et al. model
- *AFI_ODEfile.cpp*			compiled file containing ODEs for the Fogli Iseppe et al. model
- APInfo.hpp				helper file for metrics extraction
- *Drug_effects_all_\*X.dat*		parameters for simulating the drug set at different doses

### "create_datasets" folder

- *create_df_all_with_EAD_flag.py*	generates the datasets adding the info about EAD presence
- *create_training_and_test_dfs.py*	generates the training and test datasets
- *all_drugs_epi_july2020.xlsx*		file containing data related to the simulated drugs and presence of EADs in the simulations

### "RFE" folder

- *rfe_lr_bayes_optimized_repeated_n_times.py*		runs the RFE algorithm with a Logistic Regression model
- *rfe_svm_bayes_optimized_repeated_n_times.py*		runs the RFE algorithm with a Support Vector Machine model

________________________________________________________________________________________________________________

## Instructions

### List of packages required

- Numpy
- Pandas
- Multiprocessing
- Joblib
- Scikit-Learn
- Hyperopt

### Run the simulations

1. Enter in the *AFI_model* folder.

2. Create a *data* folder.

3. Run the python script named *run_AFI_model.py*.


### Create the datasets

1. Move the *all_drugs_epi_july2020.xlsx* file in the *simulations* folder.

2. The *create_df_all_with_EAD_flag.py* script expects a specific organization of the folders containing the data. In detail, the directory tree inside the *data* folder should follow the structure below:
- Sex folder (*Male* or *Female*)
- Hormone folder (*NH* for no hormones, *DHT_10* for low level of testosterone, *DHT_35* for high level of testosterone, *EF* for early follicular phase, *LF* for late follicular phase, *LU* for luteal phase)
		
Create the new folders in *data* and move the simulated data in the appropriate sex- and hormone-specific subfolder.

3. Create a copy of the *create_df_all_with_EAD_flag.py* script in each sex- and hormone-specific subfolder and run it. The script will generate a csv file containing the measured biomarkers.

4. Create a copy of the *create_training_and_test_dfs.py* script in each sex- and hormone-specific subfolder and run it. The script will split the previously created csv file in two files containing the delta biomarkers (drug value - control value) for the training (high risk and safe drugs) and test (intermediate risk) sets.


### Run the recursive feature elimination algorithm

1. Enter in the *RFE* folder.

2. Create a *data* folder.

3. Move the training csv files created during the previous step in the *data* folder.

4. Create a *results* folder.

5. Run the python script named *rfe_lr_bayes_optimized_repeated_n_times.py* if you want to train a logistic regression classifier, or *rfe_svm_bayes_optimized_repeated_n_times.py* if you want to train a support vector machine classifier.

6. Customizable variables: *n_jobs* is the maximal number of jobs (parallelization) the script will do, the range function on the same line of code (108 for lr, 109 for svm) decides how many times you want to repeat the entire RFE process for each input file, *n_iter* is the number of iterations for each optimization step (I generally used 200).

________________________________________________________________________________________________________________


## Reference

Fogli Iseppe A., Ni H., Zhu S., Zhang X., Coppini R., Yang P.C., Srivatsa U., Clancy C.E., Edwards A.G., Morotti S., Grandi E..
Sex‐specific classification of drug‐induced Torsade de Pointes susceptibility using cardiac simulations and machine learning.
Clinical Pharmacology & Therapeutics. 2021 Mar 26; TO BE COMPLETED WITH VOLUME AND PAGES ONCE PUBLISHED. doi: https://doi.org/10.1002/cpt.2240

Please, cite the above paper when using these codes.
