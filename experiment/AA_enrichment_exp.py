#!/usr/bin/python

"""AA_enrichment.py for experimental data"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import collections
from scipy import stats
import pandas as pd

from ddGs_res_exp import read_in_experiment

sys.path.append('../utlts')
from Gaussian_mixture import FitGaussianMixture


def get_sequence():
	"""collect wild-type amino acid identity for a given residue"""

	data_file = 'Nisthal_2019.xlsx' 
	df = pd.read_excel(data_file)

	d_seq = {}
	for res, wt_AA in zip(df['Position'], df['wtaa_AA_1L']):
		d_seq[res] = wt_AA
	return d_seq	

def get_aaseq_fail(residues_fail, d_seq):
	"""collect the amino acid identities of those residues that reject the null"""

	aaseq_fail = []
	for res in residues_fail:
		aaseq_fail.append(d_seq[res+1])
	return aaseq_fail



if __name__ == "__main__": 
	ddGs_exp = read_in_experiment()	#ordered by res-1
	fit_gaussian_mixture = FitGaussianMixture(ddGs_exp)

	data_file = 'Nisthal_2019.xlsx' 
	df = pd.read_excel(data_file)

	d_seq = get_sequence()

	ks = []
	pvalues = []
	for res, ddGs_res in enumerate(ddGs_exp):
		ddGs_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		if len(ddGs_19) == 17:	
			k, pvalue = stats.shapiro(ddGs_19)
			ks.append(k)
			pvalues.append(pvalue)

	count_aaseq = collections.Counter(d_seq.values())
	residues_fail = np.where(np.array(pvalues)<0.05)[0] 
	aaseq_fail = get_aaseq_fail(residues_fail, d_seq)
	count_aaseq_fail = collections.Counter(aaseq_fail)
	print count_aaseq_fail
	for aa, freq in count_aaseq_fail.most_common():
			pvalue = stats.binom.sf(freq, len(aaseq_fail), count_aaseq[aa]/float(len(d_seq)))	#not cdf
#			print aa, pvalue
			if pvalue<0.05:
				print aa, "{0:.2f}".format(freq/float(len(aaseq_fail))), "{0:.2f}".format(count_aaseq[aa]/float(len(d_seq))), "{0:.2E}".format(pvalue)
