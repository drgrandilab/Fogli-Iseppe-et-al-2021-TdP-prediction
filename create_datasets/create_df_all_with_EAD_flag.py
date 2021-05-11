# Please cite the following paper when using this model:
# Fogli Iseppe A, Ni H, Zhu S, Zhang X, Coppini R, Yang P-C, Srivatsa U,
# Clancy CE, Edwards AG, Morotti S, Grandi E (2021). Sex‐specific
# classification of drug‐induced Torsade de Pointes susceptibility using
# cardiac simulations and machine learning. Clinical Pharmacology &
# Therapeutics. doi: https://doi.org/10.1002/cpt.2240

### This script takes as input a list of folders and outputs a pandas dataframe containing several biomarkers of the cell

# import libraries

import numpy as np
import pandas as pd
import os


# define functions

def read_v_and_i(folder, drug):
    
    path = folder + '/V_and_current.%i.dat' % drug
    
    bcl = int(path.split('.')[1])
    
    v_and_i = pd.read_csv(path, header=None, sep=' ')

    v_and_i.columns = ['time','membrane_potential','Cai','Cajsr','Cansr','Nai','INa','INaL','Ito','ICaL','ICaNa','ICaK','IKr',
                   'IKs','IK1','INaCa_i','INaCa_ss','INaK','INab','IKb','IpCa','ICab','I_inj']

    last_beat = v_and_i[v_and_i.index>39900*(bcl/500)] ### use 79000 for 5 beats, 39900 for 3 beats
    
    return(last_beat)

#######################################################################################


def add_current_integrals(df, folder):
    
    ### this function calculates the integrals of the currents and adds them to the dataframe in input 

    for drug in df.index:
        data = read_v_and_i(folder, drug)
        for current in data.columns[2:-1]:
            name = current+'_integral'
            df.loc[drug, name] = abs(np.trapz(data[current], data['time']))
        df.loc[drug, 'ICaL_integral'] = df.loc[drug, 'ICaL_integral']+df.loc[drug, 'ICaNa_integral']+df.loc[drug, 'ICaK_integral']
        df.loc[drug, 'INaCa_integral'] = df.loc[drug, 'INaCa_i_integral'] + df.loc[drug, 'INaCa_ss_integral']
        df.drop(['INaCa_i_integral','INaCa_ss_integral','ICaNa_integral','ICaK_integral'], axis=1, inplace=True)

        # added min and max [Ca] for junctional and network SR
        df.loc[drug,'Cajsr_min'] = data['Cajsr'].min()
        df.loc[drug,'Cajsr_max'] = data['Cajsr'].max()
        df.loc[drug,'Cansr_min'] = data['Cansr'].min()
        df.loc[drug,'Cansr_max'] = data['Cansr'].max()

    return(df)


#########################################################################################


def read_data(folder_list, column_list=None, beat=998):
    
    beat = str(beat)
    
    if not(column_list):
        column_list = ["Stimulus_counter", "Stimulus_interval", "APD90", "APD75", "APD50", "APD30", 
                       "APD20", "Vmax", "Vmin", "dVdtmax", "plateau_potential", "INa_max",
                       "dVdt_Time_RP", "V_m60_Time_RP", "CaT_max", "CaT_min", "CaD50", "CaD80", 
                       "CaD90", "CaD_Tau", "CaD_peak_T", "dCadt_time", "qNet", "qNet_inward",
                       'Cai_integral', 'Cajsr_integral','Cansr_integral', 'Nai_integral', # added on 05/27/2020
                       'INa_integral', 'INaL_integral', 'Ito_integral', 'ICaL_integral',
                       'IKr_integral', 'IKs_integral', 'IK1_integral', 'INaK_integral',
                       'INab_integral', 'IKb_integral', 'IpCa_integral', 'ICab_integral','INaCa_integral',
                       'Cajsr_min','Cajsr_max','Cansr_min','Cansr_max'] # added on 06/01/2020
    
    cols = ['Drug','gender_type'] + column_list
    
    cell_dict = {'0':'Endo', 
                 '1':'Epi',
                 '2':'M'}

    gender_dict = {'1':'Male',
                   '2':'Female'}

    cwd = os.getcwd().split('/')
    
    ead_file = pd.read_excel('../../../all_drugs_epi_july2020.xlsx', sheet_name='EAD_'+cwd[-2].lower()+'_'+cwd[-1])

    for idx, folder in enumerate(folder_list):
        
        sub_df = pd.DataFrame(columns=cols)
        
        BCL = folder.split('.')[1]
        cell_type = folder.split('.')[3]
        gender_type = folder.split('.')[5]
        hormone_level = folder.split('.')[7]
        Cmax = folder.split('.')[9]
        drug_set = folder.split('.')[10]
        
        print(BCL,cell_type,gender_type,hormone_level,Cmax,drug_set)

        condition = BCL+'_'+Cmax

        ead_list = ead_file[condition]

        #ead_list = pd.read_csv('EAD.csv')[condition]
    
        APCaT_df = pd.read_csv(folder + '/APCaT.log.BCL.%s.CT.%s.GT.%s.HL.%s.Cmax.%s.%s.dat' 
                         % (BCL,cell_type,gender_type,hormone_level, Cmax, drug_set), 
                         sep='\s+', header=None)
        
        APCaT_df.columns = ['Population_ID','Drug_ID','Stimulus_counter498','Stimulus_interval498','APD90498',
                            'APD75498','APD50498','APD30498','APD20498','Vmax498','Vmin498','dVdtmax498',
                            'plateau_potential498','INa_max498','dVdt_Time_RP498','V_m60_Time_RP498',
                            'CaT_max498','CaT_min498','CaD50498','CaD80498','CaD90498','CaD_Tau498',
                            'CaD_peak_T498','dCadt_time498','qNet498','qNet_inward498','Stimulus_counter499',
                            'Stimulus_interval499','APD90499','APD75499','APD50499','APD30499','APD20499',
                            'Vmax499','Vmin499','dVdtmax499','plateau_potential499','INa_max499',
                            'dVdt_Time_RP499','V_m60_Time_RP499','CaT_max499','CaT_min499','CaD50499',
                            'CaD80499','CaD90499','CaD_Tau499','CaD_peak_T499','dCadt_time499','qNet499',
                            'qNet_inward499','Stimulus_counter500','Stimulus_interval500','APD90500','APD75500',
                            'APD50500','APD30500','APD20500','Vmax500','Vmin500','dVdtmax500',
                            'plateau_potential500','INa_max500','dVdt_Time_RP500','V_m60_Time_RP500',
                            'CaT_max500','CaT_min500','CaD50500','CaD80500','CaD90500','CaD_Tau500',
                            'CaD_peak_T500','dCadt_time500','qNet500','qNet_inward500','Stim_counter']
        
        
        #APCaT_sorted = APCaT_df.sort_values('Drug_ID').reset_index().drop('index', axis=1)
        APCaT_sorted = APCaT_df.set_index('Drug_ID').sort_index()
        
        APCaT_sorted = add_current_integrals(APCaT_sorted, folder)
        
        for i in APCaT_sorted.index:
            
            if ead_list[i]==0:

                values = []#["Population_ID"]
                
                for col in column_list:
                    if col[-8:]=='integral':
                        values.append(APCaT_sorted[col][i])
                    elif col=='Drug_ID':
                        values.append(APCaT_sorted[col][i])
                    elif col[:5]=='Cajsr':
                        values.append(APCaT_sorted[col][i])
                    elif col[:5]=='Cansr':
                        values.append(APCaT_sorted[col][i])
                    else:
                        values.append(APCaT_sorted[col+beat][i])

            else:
                values = [None]*45
                    
            sub_df.loc[i] = [i, gender_type]+values # 42 is a random number
    
        
        #sub_df['Population_ID'] = APCaT_df["Population_ID"].copy()

        sub_df['gender_type'] = sub_df['gender_type'].replace(gender_dict)
 
        #sub_df.sort_values('Population_ID', inplace=True)
        #sub_df.reset_index(inplace=True)
        #sub_df.drop('index', axis=1, inplace=True)
        
        #add_current_integrals(sub_df, folder) # add current integrals to dataframe    
    
        sub_df.drop(['APD20','dCadt_time','INab_integral','IKb_integral','ICab_integral'], axis=1, inplace=True)
        sub_df['AP_amplitude'] = sub_df['Vmax']-sub_df['Vmin']
        sub_df['AP_triangulation'] = sub_df['APD90']-sub_df['APD30']
    
        sub_df.drop(['Stimulus_counter','Stimulus_interval'], axis=1, inplace=True)
        sub_df.columns = list(sub_df.columns[:2]) + ['%s_BCL_%d_Cmax_%d'%(i,int(BCL),int(Cmax)) for i in sub_df.columns[2:]]
    
        if idx==0:
            df = sub_df.copy()
        else:
            sub_df.drop(['Drug', 'gender_type'], axis=1, inplace=True)
            df = pd.merge(df, sub_df, left_index=True, right_index=True)
    
    return(df)


# actual script

cwd = os.getcwd().split('/')

ans = read_data([i for i in os.listdir(os.getcwd()) if i.split('.')[0]=='BCL'], beat=498)

ans.to_csv('ead_'+cwd[-2].lower()+'_'+cwd[-1]+'.csv')

#print(os.listdir())
