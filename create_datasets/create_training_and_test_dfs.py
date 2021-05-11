### This script takes as input a pandas dataframe containing the biomarkers and outputs two dataframes of deltas, 
### one containing training drugs and one containing testing (intermediate) drugs


# import libraries

import numpy as np
import pandas as pd
import os


# define functions

def select_H_and_L_clean(df):
    
### this function gives back drugs with high and low TdP risk for which CiPA and CredibleMeds agree
    
    df_highlow = df[df['Drug'].isin([
    				'Ajmaline','Amiodarone1','Amiodarone2','Bepridil1','Bepridil2',
                 	'Bepridil_CiPA','Ceftriaxone','Cibenzoline','Cilostazol','Diazepam',
                    'Diltiazem1','Diltiazem2','Diltiazem_CiPA','Disopyramide','Dofetilide1',
                    'Dofetilide2','Dofetilide_CiPA','Donepezil','Duloxetine','Flecainide',
                    'Halofantrine','Haloperidol1','Haloperidol2','Ibutilide','Lamivudine',
                    'Linezolid','Loratadine','Methadone','Mexiletine','Mexiletine_CiPA',
                    'Mibefradil1','Mibefradil2','Mitoxantrone','Moxifloxacin','Nifedipine1',
                    'Nifedipine2','Nitrendipine1','Nitrendipine2','Pentobarbital',
                    'Phenytoin1','Phenytoin2','Prenylamine','Procainamide','Propranolol',
                    'Quinidine1','Quinidine2','Quinidine_CiPA','Ribavirin','Sitagliptin',
                    'Sotalol','Sotalol_CiPA','Sparfloxacin','Tedisamil','Telbivudine',
                    'Terodiline','Thioridazine1','Thioridazine2','Verapamil1','Verapamil2',
                    'Verapamil_CiPA'])]
    
    return(df_highlow)

############################################

def select_I_clean(df):

### this function gives back drugs with intermediate TdP risk for which CiPA and CredibleMeds agree
    
    df_inter = df[~df['Drug'].isin([
    				'Ajmaline','Amiodarone1','Amiodarone2','Bepridil1','Bepridil2',
                 	'Bepridil_CiPA','Ceftriaxone','Cibenzoline','Cilostazol','Diazepam',
                    'Diltiazem1','Diltiazem2','Diltiazem_CiPA','Disopyramide','Dofetilide1',
                    'Dofetilide2','Dofetilide_CiPA','Donepezil','Duloxetine','Flecainide',
                    'Halofantrine','Haloperidol1','Haloperidol2','Ibutilide','Lamivudine',
                    'Linezolid','Loratadine','Methadone','Mexiletine','Mexiletine_CiPA',
                    'Mibefradil1','Mibefradil2','Mitoxantrone','Moxifloxacin','Nifedipine1',
                    'Nifedipine2','Nitrendipine1','Nitrendipine2','Pentobarbital',
                    'Phenytoin1','Phenytoin2','Prenylamine','Procainamide','Propranolol',
                    'Quinidine1','Quinidine2','Quinidine_CiPA','Ribavirin','Sitagliptin',
                    'Sotalol','Sotalol_CiPA','Sparfloxacin','Tedisamil','Telbivudine',
                    'Terodiline','Thioridazine1','Thioridazine2','Verapamil1','Verapamil2',
                    'Verapamil_CiPA'])]
    
    return(df_inter)

############################################

def fill_na_ead(df):

### this function takes as input a dataframe with NaNs and return the dataframe filled with the
### maximum of minimum of the TdP+ group value for each column, depending on the fact that the
### average value of TdP+ group is larger or smaller than the average value of TdP- group,
### respectively.

	p = df.groupby('TdP').mean()

	for col in df.columns[2:-18]:

		#print(col)

		if (p.loc[1,col]>p.loc[0,col]):
			df[col].fillna(df.groupby('TdP').max().loc[1,col], inplace=True)

		else:
			df[col].fillna(df.groupby('TdP').min().loc[1,col], inplace=True)

	return(df)

############################################

### ACTUAL SCRIPT

file = [i for i in os.listdir() if i[-3:]=='csv'][0]

data = pd.read_csv(file, index_col=0)
# Haibo modified the folder to locate ''../../../all_drugs_epi_july2020.xlsx''
# IC50_drugs = pd.read_excel('/media/afogliiseppe/DATADRIVE0/TdP_project_2020/july2020_model/ALL_EPI/all_drugs_epi_july2020.xlsx', nrows=99)
IC50_drugs = pd.read_excel('../../../all_drugs_epi_july2020.xlsx', nrows=99)

drugs_dict = {k:v for k,v in IC50_drugs['Drug name'].to_dict().items()}

data['Drug'] = data['Drug'].replace(drugs_dict)
# Haibo modified the folder to locate ''../../../all_drugs_epi_july2020.xlsx''
# ead_excel = pd.read_excel('/media/afogliiseppe/DATADRIVE0/TdP_project_2020/july2020_model/ALL_EPI/all_drugs_epi_july2020.xlsx', sheet_name='EAD'+file[3:-4])
ead_excel = pd.read_excel('../../../all_drugs_epi_july2020.xlsx', sheet_name='EAD'+file[3:-4])

cols = ['EAD_'+name for name in ead_excel.iloc[:,2:].columns]

data[cols] = ead_excel.iloc[:,2:]

data[['IC50_hERG','IC50_ICaL','IC50_INa']] = IC50_drugs[['IC50 hERG','IC50 IcaL','IC50 Ina']].copy()

# using IC50s over free concentration AFI 07/30/2020
data['IC50_hERG'] /= IC50_drugs['Free Cmax (nM)']
data['IC50_ICaL'] /= IC50_drugs['Free Cmax (nM)']
data['IC50_INa'] /= IC50_drugs['Free Cmax (nM)']

data['ratio'] = data['IC50_ICaL']/data['IC50_hERG']
data['TdP'] = IC50_drugs['TdP_risk'].copy()

for initi in ['^dVdt_','^V_m60','^CaD90','^CaD50','^CaD_peak','qNet','Ca.sr_']:
    data.drop(list(data.filter(regex=initi).columns), axis=1, inplace=True)

data_filled = fill_na_ead(data)

data_filled.iloc[:, 2:-18] = data_filled.iloc[:, 2:-18] - data_filled.iloc[0, 2:-18]

data_filled.drop(0, inplace=True)

### NOTE: I am going to save only the dataframe with condition-specific EAD flags because it was the best performing in
###       the previous feature selection (no difference with a generic EAD flag because not selected anyway)

training = select_H_and_L_clean(data_filled).drop('EAD_Any', axis=1).dropna() # TdP+ and TdP- drugs (safe and high risk)

intermediate = select_I_clean(data_filled).drop(['EAD_Any','TdP'], axis=1).dropna() # Drugs with intermediate risk

training.to_csv(file[:-4]+'_training.csv')
intermediate.to_csv(file[:-4]+'_intermediate.csv')