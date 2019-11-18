#!/usr/bin/python

help_msg = 'fit GaussianMixture of any number of Gaussians using same binning as Tokuriki et al. (2007)'

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy import stats


aminoacids = ['GLY', 'ALA', 'LEU', 'VAL', 'ILE', 'PRO', 'ARG', 'THR', 'SER', 'CYS', 'MET', 'LYS', 'GLU', 'GLN', 'ASP', 'ASN', 'TRP', 'TYR', 'PHE', 'HIS']

class FitGaussianMixture():
	#same conditions as Tokuriki 2007
	binwidth = 1 #kcal/mol
	min_ddG = -10
	max_ddG = 15 

	def __init__(self, ddGs):
		self.hist, self.bins = np.histogram(self.preprocess(ddGs), bins = np.arange(self.min_ddG, self.max_ddG + self.binwidth, self.binwidth), density=True)
		dist = (self.bins[1]-self.bins[0])/2.
		self.x = self.bins[:-1]+dist		 

	def preprocess(self, ddGs):
		ddGs = ddGs[~np.isnan(ddGs)]
		ddGs[ddGs < self.min_ddG] = self.min_ddG
		ddGs[ddGs > self.max_ddG] = self.max_ddG
		return ddGs

	def parse_params(self, params):
		weights = params[:self.nGaussians-1]	#make sure last weight is normalized!
		means = params[self.nGaussians - 1:self.nGaussians*2 - 1]
		stds = params[self.nGaussians*2 - 1:self.nGaussians*3 - 1]
		return weights, means, stds

	def MLE_func(self, params):	#can't pass x for minimize function
		weights, means, stds = self.parse_params(params)
		LLstd = params[-1]
	
		yPred = 0
		for weight, mean, std in zip(weights, means, stds):
			yPred += weight*stats.norm.pdf(self.x, loc=mean, scale=std)
		yPred += (1-np.sum(weights))*stats.norm.pdf(self.x, loc=means[-1], scale=stds[-1])
		negLL = -np.sum( stats.norm.logpdf(self.hist, loc=yPred, scale=LLstd ) )
		return negLL 
	
	def MLE(self, coeffs):
		p0 = list(coeffs)	#lower LL
		p0.append(0.02)	#LLstd
		results = minimize(self.MLE_func, p0, method='Nelder-Mead', options={'maxfev': 100000, 'disp':True})	#Nelder-Mead can't do bounds
		return results.fun, results.x

	def func2(self, coeffs):
		nGaussian = (len(coeffs)+1)/3
		weights = coeffs[:nGaussian-1]	#make sure last weight is normalized!
		means = coeffs[nGaussian - 1:nGaussian*2 - 1]
		stds = coeffs[nGaussian*2 - 1:nGaussian*3 - 1]
		LLstd = coeffs[-1]

		yPred = 0
		for weight, mean, std in zip(weights, means, stds):
			yPred += weight*stats.norm.pdf(self.x, loc=mean, scale=std)
		yPred += (1-np.sum(weights))*stats.norm.pdf(self.x, loc=means[-1], scale=stds[-1])

		negLL = -np.sum( stats.norm.logpdf(self.hist, loc=yPred, scale=LLstd ) )	#this sd is diff!

		return negLL 

	def prep_run(self):
		self.p0 = np.zeros(self.nGaussians*3 - 1)
		self.p0.fill(1)
		self.p0[:self.nGaussians-1] = 1/float(self.nGaussians)
	
		bounds = []
		for i in range(len(self.p0)):
			if i < self.nGaussians-1:
				bounds.append((0, 1))
			elif i >= self.nGaussians-1 and i < self.nGaussians*2 - 1:
				bounds.append((-np.inf, np.inf))
			else:
				bounds.append((0, np.inf))
		self.bounds = zip(*bounds)

	def run(self, nGaussians):
		self.nGaussians = nGaussians
		self.prep_run()

		LSE_coeffs, matcov = curve_fit(func, self.x, self.hist, self.p0, bounds=self.bounds)	#least-squares approach	#need bounds!
		negLL, MLE_coeffs = self.MLE(LSE_coeffs)
#		print LSE_coeffs
#		print MLE_coeffs
		return negLL, MLE_coeffs

	def special_run(self, coeffs):
		results2 = minimize(self.func2, coeffs, method='Nelder-Mead', options={'maxfev': 100000, 'disp':True})
#		print results.x

def func(x, *params):	#need this cuz curve_fit has specific syntax
	nGaussian = (len(params)+1)/3
	weights = params[:nGaussian-1]	#make sure last weight is normalized!
	means = params[nGaussian - 1:nGaussian*2 - 1]
	stds = params[nGaussian*2 - 1:]

	dist = 0
	for weight, mean, std in zip(weights, means, stds):
		dist += weight*stats.norm.pdf(x, loc=mean, scale=std)
	dist += (1-np.sum(weights))*stats.norm.pdf(x, loc=means[-1], scale=stds[-1])

	return dist 

def get_protein_ddG(df_protein):
	df_protein = df_protein[aminoacids]
	ddGs = df_protein.values
	ddGs = ddGs.astype(float)
	return ddGs

