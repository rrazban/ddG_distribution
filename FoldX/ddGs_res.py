#!/usr/bin/python

help_msg = 'determine whether ddGs per residue can be fit to a Gaussian (Figure 1 in paper)'

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

sys.path.append('../')
from Gaussian_mixture import FitGaussianMixture, get_protein_ddG


proteins = ['PON', 'lipase', '1BTL', 'CAII', 'DHFR', 'Rnase H', 'Myoglobin', 'Snase', 'human Lysozyme', 'lysozyme', 'Rnase', 'Barnase', 'AcP', 'Ubiqitin', 'Protein G', 'Cro_repressor']	#input in data file
pretty_proteins = ['PON', 'Lip', 'TEM1', 'CAII', 'DHFR', 'RNaH', 'Myo', 'SNase', 'huLys', 'henLys', 'RNaA', 'Barn', 'AcP', 'Ubiq', 'proG', 'Cro']	#output in my figure


def goodness_of_fit(ddGs):
	fit_gaussian_mixture = FitGaussianMixture(ddGs)
	ks = []
	pvalues = []
	for res, ddGs_res in enumerate(ddGs):
		ddGs_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		k, pvalue = stats.shapiro(ddGs_19)

		ks.append(k)
		pvalues.append(pvalue)
	return ks, pvalues
	
def count_fails(pvalues, num_ddGs):
	low = len(np.where(pvalues<0.05)[0])/float(len(pvalues))*100
	high = len(np.where(pvalues<0.001)[0])/float(len(pvalues))*100
	return low, high

def plotout(x, bar1, bar2):
	plt.figure(figsize=(10, 4.5))
	plt.bar(x-0.1, bar1, width=0.2, label='p-value < 0.05')
	plt.bar(x+0.1, bar2, width=0.2, label='p-value < 0.001')
	plt.xticks(x, pretty_proteins)# rotation = 20)
	plt.ylabel('% of residues significantly different from Gaussian')
	plt.legend()
	plt.show()


if __name__ == '__main__':
	data_file = 'Tokuriki_2007.xlsx'
	df = pd.read_excel(data_file)

	bar1s = []
	bar2s = []
	for p, protein_name in enumerate(proteins):
		protein = df[df['Name']==protein_name]
		ddGs = get_protein_ddG(protein)

		ks, pvalues = goodness_of_fit(ddGs)
		bar1, bar2 = count_fails(np.array(pvalues), len(ddGs))
		bar1s.append(bar1)
		bar2s.append(bar2)

	plotout(np.arange(len(proteins)), bar1s, bar2s)
