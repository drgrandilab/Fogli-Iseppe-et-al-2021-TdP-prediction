# -*- coding: utf-8 -*-
# Please cite the following paper when using this model:
# Fogli Iseppe A, Ni H, Zhu S, Zhang X, Coppini R, Yang P-C, Srivatsa U,
# Clancy CE, Edwards AG, Morotti S, Grandi E (2021). Sex‐specific
# classification of drug‐induced Torsade de Pointes susceptibility using
# cardiac simulations and machine learning. Clinical Pharmacology &
# Therapeutics. doi: https://doi.org/10.1002/cpt.2240

# NOTE: you need to create an empty "data" folder before running the code

from subprocess import call
import math
import sys
import os
from joblib import Parallel, delayed
import multiprocessing
import numpy as np


def run_Normal(celltype, typegender,Homlevel, PopulID,  PopulPara, drugID, drugeffect, BCL,  f):
	P=['./AFI_ODEfile.out']
	P.append(str(celltype))
	P.append(str(type_gender))
	P.append(str(Homlevel))
	P.append(str(PopulID))


	for i in range(len(PopulPara)):
		P.append(str(PopulPara[i]))
	for i in range(len(drugeffect)):
		P.append(str(drugeffect[i]))
	P.append(str(drugID))
	P.append(str(BCL))
	call(P, stdout=f)


def simulate_drugs(i):
	print(i)
	drugID = i;
	drugeffect = drug_effect_list[i]

	run_Normal(cellType, type_gender, Homlevel, PopulID,  PopulPara, drugID, drugeffect, BCL,  files)
	call('mv V_and_current.%d.dat '%drugID + path + name_folder, shell=True)


path = os.getcwd()

num_cores = multiprocessing.cpu_count() - 1

PopulPara = [1]*13; # parameters to generate a population of models
cellType = 1;
PopulID = 0;

for BCL in [500.0, 1000.0, 2000.0]: # basic cycle length
	for Cmax in [1, 2, 3, 4]: # Drug concentration

		type_gender = 1

		# 0 for baseline, 1 for DHT 10 nM, 2 for DHT 35 nM
		for Homlevel in [0]:

			drug_effect_list = np.loadtxt('Drug_effects_all_%dX.dat'%Cmax)

			name_folder = '/data/BCL.%d.CT.%d.GT.%d.HL.%d.Cmax.%d.ALL'%(BCL,cellType,type_gender,Homlevel,Cmax)
			os.mkdir(path + name_folder)

			files=open("APCaT.log.BCL.%d.CT.%d.GT.%d.HL.%d.Cmax.%d.ALL.dat"%(BCL,cellType,type_gender,Homlevel,Cmax), "w+")

			Parallel(n_jobs=num_cores, prefer='threads')(delayed(simulate_drugs)(i) for i in range(len(drug_effect_list)))

			call('mv *.txt summary.dat APCaT.log.BCL.%d.CT.%d.GT.%d.HL.%d.Cmax.%d.ALL.dat '%(BCL,cellType,type_gender,Homlevel,Cmax)+ path + name_folder, shell=True)


		type_gender = 2

		# 0 for baseline, 1 for early follicular phase, 2 for late follicular phase, 3 for luteal phase
		for Homlevel in [0]:

			drug_effect_list = np.loadtxt('Drug_effects_all_%dX.dat'%Cmax)

			name_folder = '/data/BCL.%d.CT.%d.GT.%d.HL.%d.Cmax.%d.ALL'%(BCL,cellType,type_gender,Homlevel,Cmax)
			os.mkdir(path + name_folder)

			files=open("APCaT.log.BCL.%d.CT.%d.GT.%d.HL.%d.Cmax.%d.ALL.dat"%(BCL,cellType,type_gender,Homlevel,Cmax), "w+")

			Parallel(n_jobs=num_cores, prefer='threads')(delayed(simulate_drugs)(i) for i in range(len(drug_effect_list)))

			call('mv *.txt summary.dat APCaT.log.BCL.%d.CT.%d.GT.%d.HL.%d.Cmax.%d.ALL.dat '%(BCL,cellType,type_gender,Homlevel,Cmax)+ path + name_folder, shell=True)
