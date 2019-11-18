#!/usr/bin/python

help_msg = 'compare overall ddG distribution with N Gaussian model and bi-Gaussian fit (Figure 2 in paper)'

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from ddGs_res import proteins, pretty_proteins

sys.path.append('../')
from Gaussian_mixture import FitGaussianMixture, get_protein_ddG


def generate_dist(ddGs, fit_gaussian_mixture):
	yPred = 0
	for res, ddGs_res in enumerate(ddGs):	#not efficient but want to make sure preprocess done correctly
		ddG_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		fit_mu, fit_std=stats.norm.fit(ddG_19)
		yPred += stats.norm.pdf(fit_gaussian_mixture.x, loc=fit_mu, scale=fit_std)	

	yPred /= float(len(ddGs))	#normalize
	return yPred

def KS_pvalue_bound(pvalue, num1, num2):	#not consistent with how python calculates pvalue, seems to be distr. specific
	c = np.sqrt(-0.5*np.log(pvalue))
	k = c*np.sqrt((num1+num2)/float(num1*num2))
	return k

def KS_pvalue(k, num1, num2):
	return np.exp(-2*k**2 * num1*num2/(num1 + num2))

def run_gaussian_mixture(fit_gaussian_mixture):
	nGaussian = 2
	negLL, MLE_coeffs = fit_gaussian_mixture.run(nGaussian)
	weights, means, stds = fit_gaussian_mixture.parse_params(MLE_coeffs)

	Gaussian = 0
	x = fit_gaussian_mixture.x
	for weight, mean, std in zip(weights, means, stds):
		Gaussian += weight*stats.norm.pdf(x, loc=mean, scale=std)
	Gaussian += (1-np.sum(weights))*stats.norm.pdf(x, loc=means[-1], scale=stds[-1])

	k2, pvalue2 = stats.ks_2samp(fit_gaussian_mixture.hist, Gaussian)
	return k2, pvalue2, Gaussian

def plotout(x, bar1, bar2):
	plt.figure(figsize=(10, 4.5))
	plt.bar(x-0.1, bar1, width = 0.2, label='N-Gaussian model')
	plt.bar(x+0.1, bar2, width = 0.2, label='bi-Gaussian fit')
	plt.axhline(KS_pvalue_bound(0.05, len(x), len(x)), label = 'p-value = 0.05', linestyle = ':', color = 'k')
	plt.axhline(KS_pvalue_bound(0.001, len(x), len(x)), label = 'p-value = 0.001', linestyle = '--', color = 'k')

	plt.xticks(x, pretty_proteins)# rotation = 20)
	plt.ylabel('K-S test statistic')
	plt.legend()
	plt.show()

if __name__ == '__main__':
	data_file = 'Tokuriki_2007.xlsx'
	df = pd.read_excel(data_file)

	ks = []
	ks2 = []
	for protein_name in proteins:
		protein = df[df['Name']==protein_name]

		ddGs = get_protein_ddG(protein)
		fit_gaussian_mixture = FitGaussianMixture(ddGs)
		ddGs_overall = fit_gaussian_mixture.preprocess(ddGs)	#destroys shape

		yPred = generate_dist(ddGs, fit_gaussian_mixture)

		k, pvalue = stats.ks_2samp(fit_gaussian_mixture.hist, yPred)

		ks.append(k)
		ks2.append(run_gaussian_mixture(fit_gaussian_mixture)[0])

	plotout(np.arange(len(proteins)), ks, ks2)



