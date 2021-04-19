Fogli Iseppe et al. male and female human epicardial ventricular cardiomyocyte models.

The parent model developed by Yang and Clancy (Yang, P.-C. & Clancy, C. E. In silico Prediction of Sex-Based Differences in Human Susceptibility to Cardiac Ventricular Tachyarrhythmias. Front. Physiol. 3, 1–12 (2012)) was modified by:
	- increasing the maximal transport rate of the Na+-Ca2+ exchanger (NCX) by 15% in the female model
	- removing the original female-to-male difference in SERCA formulation
	- removing the original female-to-male difference in Na+/K+-ATPase formulation
The updated model better recapitulates observed functional sex differences in Ca2+ handling.

_______________________________________________________________________________________________________________

Contents:

readme.txt		 	this file

"AFI_model" folder

	run_AFI_model.py		runs the simulations with the baseline male and female models
	AFI_ODEfile.cpp			source file containing ODEs for the Fogli Iseppe et al. model
	AFI_ODEfile.cpp			compiled file containing ODEs for the Fogli Iseppe et al. model
	APInfo.hpp			helper file for metrics extraction
	Drug_effects_all_*X.dat 	parameters for simulating the drug set at different doses

________________________________________________________________________________________________________________


Reference:
Fogli Iseppe A., Ni H., Zhu S., Zhang X., Coppini R., Yang P.C., Srivatsa U., Clancy C.E., Edwards A.G., Morotti S., Grandi E..
Sex‐specific classification of drug‐induced Torsade de Pointes susceptibility using cardiac simulations and machine learning.
Clinical Pharmacology & Therapeutics. 2021 Mar 26; TO BE COMPLETED WITH VOLUME AND PAGES ONCE PUBLISHED. doi: https://doi.org/10.1002/cpt.2240

Please, cite the above paper when using these codes.
