#!/usr/bin/python

help_msg = 'ddGs_res.py for one protein (Figure S1 in paper)'

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

sys.path.append('../')
from Gaussian_mixture import FitGaussianMixture, get_protein_ddG


def plotout_res(ddGs, k, pvalue, x, gaussian, bins):
	plt.hist(ddGs, bins, density = True)	
	plt.plot(x, gaussian, color='k', label = 'W = {0:.2f}\np-value = {1:.2E}'.format(k, pvalue))
	plt.title('first residue of dihydrofolate reductase')
	plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
	plt.ylabel('probability')
	plt.legend()
	plt.show()
#	sys.exit()

def plotout(ks):
	plt.hist(ks)
	plt.xlabel('K-S statistic')
	plt.title('all DHFR residues')
	plt.ylabel('frequency')
	plt.legend()
	plt.show()


if __name__ == '__main__':
	data_file = 'Tokuriki_2007.xlsx'
	df = pd.read_excel(data_file)

	protein = df[df['Name']=='DHFR']

	ddGs = get_protein_ddG(protein)
	length = len(ddGs)
	fit_gaussian_mixture = FitGaussianMixture(ddGs)

	x = fit_gaussian_mixture.x
	ASAs = (protein['ASA']).values

	for res, ddGs_res in enumerate(ddGs):
		ddGs_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		fit_mu, fit_std=stats.norm.fit(ddGs_19)
		k, pvalue = stats.shapiro(ddGs_19)

		if True:
			plotout_res(ddGs_19, k, pvalue, x, stats.norm.pdf(x, fit_mu, fit_std), fit_gaussian_mixture.bins)	
			sys.exit()
