#!/usr/bin/python

"""ddGs_res.py for experimental data (Figure 3)"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from Bio.SeqUtils import seq1

sys.path.append('../utlts')
from Gaussian_mixture import FitGaussianMixture, aminoacids

sys.path.append('../FoldX')
from ddGs_res_DHFR import plotout_res


def parse(res):
	"""return residue number and amino acid identity"""	

	if len(res)==3:
		return int(res[1]), res[2]
	elif len(res)==4:
		return int(res[1:3]), res[3]

def read_in_experiment():
	""" 
	Read in raw data from Nisthal_2019.xlsx, making sure
	it matches FoldX read in
  
    Returns: 
	array-like   
		organized ddG data by residue and amino acid type

    """

	AAs = [seq1(aa_3) for aa_3 in aminoacids]

	data_file = 'Nisthal_2019.xlsx'	
	df2 = pd.read_excel(data_file)
	ddGs = np.zeros((56, 20))
	ddGs.fill(np.nan)	#match Tokuriki dataset read in
	for res, ddG in zip(df2['MUT_LBL'], df2['ddG(mAvg)_mean']):
		res_num, AA = parse(res)
		try:
			ddG = float(ddG)
			if ddG!=-4:	#Nisthal includes the black squares in fig 2 as -4 kcal/mol
				ddGs[res_num-1][AAs.index(AA)] = -float(ddG)
		except:
			pass
	return ddGs

def get_bars(all_pvalues, some_pvalues):
	"""obtain number of residues that reject the null"""

	bars = []
	for threshold in [0.05, 0.001]:
		bar1 = []
		for pvalues in [all_pvalues, some_pvalues]:
			bar1.append(len(np.where(pvalues < threshold)[0]))
#			print list(np.where(pvalues < threshold)[0])
		bars.append(bar1)
	return bars

def plotout(x, bar1, bar2):
	plt.title('$\\beta$1 domain of protein G')
	plt.bar(x-0.1, bar1, width=0.2, label='p-value < 0.05')
	plt.bar(x+0.1, bar2, width=0.2, label='p-value < 0.001', hatch='//')
	plt.xticks(x, ['>2', '17 (max)'])# rotation = 20)
	plt.xlabel('number of measurements per residue')
	plt.ylabel('# of residues significantly different from Gaussian')
	plt.legend()
	plt.show()



if __name__ == '__main__':
	ddGs_exp = read_in_experiment()
	fit_gaussian_mixture = FitGaussianMixture(ddGs_exp)
	ddGs_all = fit_gaussian_mixture.preprocess(ddGs_exp)
#	print len(ddGs_all), len(ddGs_all)/(56.*19)
	x = fit_gaussian_mixture.x

	ks = []
	pvalues = []
	num_mutants = []
	for res, ddGs_res in enumerate(ddGs_exp):
		ddGs_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		if len(ddGs_19) > 2:	#minimal number of datapoints to get a variance
			num_mutants.append(len(ddGs_19))
			fit_mu, fit_std=stats.norm.fit(ddGs_19)
			k, pvalue = stats.shapiro(ddGs_19)
			ks.append(k)
			pvalues.append(pvalue)

			if False:
				plotout_res(ddGs_19, k, pvalue, x, stats.norm.pdf(x, fit_mu, fit_std), fit_gaussian_mixture.bins)	
				sys.exit()
	pvalues = np.array(pvalues)
	num_mutants = np.array(num_mutants)
	bars = get_bars(pvalues, pvalues[np.where(num_mutants==17)[0]])
	plotout(np.array([0,1]), bars[0], bars[1]) 
