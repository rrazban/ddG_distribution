#!/usr/bin/python

help_msg = 'analyze amino acids of ddGs_res.py residues that reject the null, compare to overall proteins amino acid usage (Figure S2 in paper)'

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import collections
from scipy import stats
import pandas as pd

sys.path.append('../')
from Gaussian_mixture import get_protein_ddG
from ddGs_res import goodness_of_fit, proteins, pretty_proteins


def get_aaseq_fail(residues_fail, protein):
	residues_ref = list(protein['No'])

	aaseq_fail = []
	for res in residues_fail:
		res_ref = residues_ref[res]
		row = protein[protein['No']==res_ref]
		aaseq_fail.extend(list(row['A.A']))
	return aaseq_fail

def prepare_output(AAs, count_enriched_AAs):
	d_translate = {u'GLY': 'G', u'LYS': 'K', u'CYS': 'C', u'ASP': 'D', u'GLN': 'Q', u'SER': 'S', u'PRO': 'P', u'ALA': 'A', u'HIS': 'H', u'ARG': 'R', u'TYR': 'Y', u'VAL': 'V', u'MET': 'M', u'ILE': 'I', u'PHE': 'F', u'TRP': 'W', u'THR': 'T', u'ASN': 'N', u'LEU': 'L', u'GLU': 'E'}
	inverted_dict = dict([[v,k] for k,v in d_translate.items()])

	output = []
	for AA in AAs:
		AA_3 = inverted_dict[AA]
		if AA_3 in count_enriched_AAs:
			output.append(count_enriched_AAs[AA_3])
		else:
			output.append(0)
	return output

def plotout(xlabels, output):
	xs = range(len(output))
	xs.extend([4, 12, 17])
	xlabels.extend(['\nhydrophobic', '\npolar', '\ncharged']) 
	plt.bar(range(len(output)), output)
	plt.xticks(xs, xlabels) 
	plt.ylabel('# of proteins with enriched amino acid')
	plt.show()


if __name__ == "__main__": 
	data_file = 'Tokuriki_2007.xlsx'
	df = pd.read_excel(data_file)

	AAs = 'AVILMFYWCGPSTNQRHKDE'	#ordered

	enriched_AAs = []
	for p, protein_name in enumerate(proteins):
		protein = df[df['Name']==protein_name]
		aaseq = list(protein['A.A'])
		count_aaseq = collections.Counter(aaseq)
		ddGs = get_protein_ddG(protein)

		ks, pvalues = goodness_of_fit(ddGs)
		residues_fail = np.where(np.array(pvalues)<0.05)[0]

		aaseq_fail = get_aaseq_fail(residues_fail, protein)	
		count_aaseq_fail = collections.Counter(aaseq_fail)
#		print pretty_proteins[p], count_aaseq

		for aa, freq in count_aaseq_fail.most_common():
			pvalue = stats.binom.sf(freq, len(aaseq_fail), count_aaseq[aa]/float(len(aaseq)))	#not cdf
			if pvalue<0.05:
#				print aa, "{0:.2f}".format(freq/float(len(aaseq_fail))), "{0:.2f}".format(count_aaseq[aa]/float(len(aaseq))), "{0:.2E}".format(pvalue)
				enriched_AAs.append(aa)
#		print pretty_proteins[p], count_aaseq_fail

	count_enriched_AAs = collections.Counter(enriched_AAs)
#	print count_enriched_AAs, len(count_enriched_AAs)

	output = prepare_output(AAs, count_enriched_AAs)
	plotout(list(AAs), output)
