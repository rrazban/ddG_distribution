#!/usr/bin/python

"""contact_dist.py for experimental data"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from ddGs_res_exp import read_in_experiment

sys.path.append('../utlts')
from Gaussian_mixture import FitGaussianMixture

sys.path.append('../FoldX')
from contact_dist import get_contacts


def plotout(protein_name, rejected_contacts, all_contacts):
	plt.title(protein_name)
	plt.boxplot([rejected_contacts, all_contacts], showmeans=True)#, label = 'all')
	plt.xlabel('residues')
	plt.xticks([1, 2], ['reject the null', 'all'])
	plt.ylabel('frequency of residues')
	k, pvalue = stats.ks_2samp(all_contacts, rejected_contacts)	#only CAII significant
	plt.legend(title = 'KS = {0:.2f} ({1:.2E})'.format(k, pvalue))
	plt.show()	


if __name__ == "__main__": 
	ddGs_exp = read_in_experiment()	#ordered by res-1
	fit_gaussian_mixture = FitGaussianMixture(ddGs_exp)

	ks = []
	pvalues = []
	residues_ref = []
	for res, ddGs_res in enumerate(ddGs_exp):
		ddGs_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		if len(ddGs_19) == 17:	
			k, pvalue = stats.shapiro(ddGs_19)
			ks.append(k)
			pvalues.append(pvalue)
			residues_ref.append(res+1)	#account for res-1

	residues_fail = np.where(np.array(pvalues)<0.05)[0] 

	d_contacts = get_contacts('1pga')
	all_contacts = d_contacts.values()
	rejected_contacts = []	
	for res in residues_fail:
		res_ref = residues_ref[res]
		rejected_contacts.append(d_contacts[res_ref])
	plotout('GB1', rejected_contacts, all_contacts)
