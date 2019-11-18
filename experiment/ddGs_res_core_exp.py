#!/usr/bin/python

help_msg = 'ddGs_res_core.py for experimental data'

import sys
import numpy as np
from scipy import stats

from ddGs_res_exp import read_in_experiment

sys.path.append('../')
from Gaussian_mixture import FitGaussianMixture
sys.path.append('../FoldX')
from ddGs_res_core import plotout, count_fails


def goodness_of_fit(ddGs):
	fit_gaussian_mixture = FitGaussianMixture(ddGs)
	ks = []
	pvalues = []
	for res, ddGs_res in enumerate(ddGs):
		ddGs_19 = fit_gaussian_mixture.preprocess(ddGs_res)
		if len(ddGs_19)==17:
			k, pvalue = stats.shapiro(ddGs_19)
			ks.append(k)
			pvalues.append(pvalue)
	return ks, pvalues
	
def transform(d_ASAs):
	ASAs = np.zeros(len(d_ASAs))
	d_surface = {'core': 0, 'surface':1}
	for res, loc in d_ASAs.iteritems():
		ASAs[res-1] = d_surface[loc]
	return ASAs

if __name__ == '__main__':
	ddGs_exp = read_in_experiment()
	d_ASAs = {1: 'surface', 2: 'surface', 3: 'core', 4: 'surface', 5: 'core', 6: 'surface', 7: 'core', 8: 'surface', 9: 'core', 10: 'surface', 11: 'surface', 12: 'surface', 13: 'surface', 14: 'surface', 15: 'surface', 16: 'surface', 17: 'surface', 18: 'core', 19: 'surface', 20: 'core', 21: 'surface', 22: 'surface', 23: 'surface', 24: 'surface', 25: 'surface', 26: 'core', 27: 'surface', 28: 'surface', 29: 'surface', 30: 'core', 31: 'surface', 32: 'surface', 33: 'surface', 34: 'core', 35: 'surface', 36: 'surface', 37: 'surface', 38: 'surface', 39: 'core', 40: 'surface', 41: 'surface', 42: 'surface', 43: 'surface', 44: 'surface', 45: 'surface', 46: 'surface', 47: 'surface', 48: 'surface', 49: 'surface', 50: 'surface', 51: 'core', 52: 'core', 53: 'core', 54: 'core', 55: 'surface', 56: 'surface'}	#taken from my analysis on 1pga; same ASA thresholds used as in Tokuriki 2007
	ASAs = transform(d_ASAs)

	bar1s = []
	bar2s = []

	ks, pvalues = goodness_of_fit(ddGs_exp)
	bar1, bar2 = count_fails(np.array(pvalues), ASAs)
	print bar1, bar2

	plotout(np.array([1]), bar1, bar2)

