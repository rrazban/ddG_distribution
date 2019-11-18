#!/usr/bin/python

help_msg = 'parse ddGs_res.py into core residues, compare to %core residues in protein (Figure S2 in paper)'

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from ddGs_res import proteins, pretty_proteins, goodness_of_fit

sys.path.append('../')
from Gaussian_mixture import FitGaussianMixture, get_protein_ddG

def count_fails(pvalues, ASAs):
	fails_res = np.where(pvalues<0.05)[0]
	core_fail = 0
	for res in fails_res:
		if ASAs[res]<0.25:
			core_fail+=1
	
	ref_core = len(np.where(ASAs<0.25)[0])/float(len(ASAs)) * 100
	try:
		return core_fail/float(len(fails_res)) * 100, ref_core
	except:
		return 0, ref_core

def plotout(x, bar1, bar2):
	plt.figure(figsize=(10, 4.5))

	plt.bar(x-0.1, bar1, width=0.2, label='p-value < 0.05, core')
	plt.bar(x+0.1, bar2, width=0.2, label='protein core')
	plt.xticks(x, pretty_proteins)# rotation = 20)
	plt.ylabel('% of residues')
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
		ASAs = (protein['ASA']).values

		ks, pvalues = goodness_of_fit(ddGs)
		bar1, bar2 = count_fails(np.array(pvalues), ASAs)
		bar1s.append(bar1)
		bar2s.append(bar2)

	plotout(np.arange(len(proteins)), bar1s, bar2s)

