#!/usr/bin/python

"""ddGs.py for experimental data (Figure S4 in paper)"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from ddGs_res_exp import read_in_experiment

sys.path.append('../utlts')
from Gaussian_mixture import FitGaussianMixture

sys.path.append('../FoldX')
from ddGs_res import proteins, pretty_proteins
from ddGs import run_gaussian_mixture, generate_dist, KS_pvalue


def parse(ddGs):
	"""identify those residues with less than the maximum number (17) of mutant ddGs measured"""

	del_list = []
	for res, ddG_res in enumerate(ddGs):
		ddGs_17 = fit_gaussian_mixture.preprocess(ddG_res)
		if len(ddGs_17)!=17:
			del_list.append(res)
 	return del_list

def plotout(ddGs_exp, fit_gaussian_mixture):
	plt.title('$\\beta$1 domain of protein G')
	plt.hist(ddGs_exp[~np.isnan(ddGs_exp)], fit_gaussian_mixture.bins, density = True)	
	plt.plot(fit_gaussian_mixture.x, yPred, color='k', label = 'N-Gaussian\nD = {0:.2f}\np-value = {1:.2E}'.format(k, KS_pvalue(k, len(fit_gaussian_mixture.x), len(fit_gaussian_mixture.x))))
	plt.plot(fit_gaussian_mixture.x, gaussian, color='r', label = 'bi-Gaussian\nD = {0:.2f}\np-value = {1:.2E}'.format(ks2, KS_pvalue(ks2, len(fit_gaussian_mixture.x), len(fit_gaussian_mixture.x))))
	plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
	plt.ylabel('probability')
	plt.legend()
	plt.show()

if '__main__':
	ddGs_exp = read_in_experiment()
		
	fit_gaussian_mixture = FitGaussianMixture(ddGs_exp)
	del_list = parse(ddGs_exp)
	ddGs_exp = np.delete(ddGs_exp, del_list, axis = 0)

	yPred = generate_dist(ddGs_exp, fit_gaussian_mixture)

	k, pvalue = stats.ks_2samp(fit_gaussian_mixture.hist, yPred)
	ks2, pvalue2, gaussian = run_gaussian_mixture(fit_gaussian_mixture)

	plotout(ddGs_exp, fit_gaussian_mixture)
