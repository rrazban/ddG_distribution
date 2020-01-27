#!/usr/bin/python

"""analyze amino acids of ddGs_res.py residues that reject the null, compare to overall proteins amino acid usage (Figure S2 in paper)"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import collections
from scipy import stats
import pandas as pd
from Bio.SeqUtils import seq3

from ddGs_res import goodness_of_fit, proteins, pretty_proteins

sys.path.append('../utlts')
from Gaussian_mixture import get_protein_ddG


def get_aaseq_fail(residues_fail, protein):
	"""collect the amino acid identities of those residues that reject the null"""

	residues_ref = list(protein['No'])

	aaseq_fail = []
	for res in residues_fail:
		res_ref = residues_ref[res]
		row = protein[protein['No']==res_ref]
		aaseq_fail.extend(list(row['A.A']))
	return aaseq_fail

def prepare_output(AAs, count_enriched_AAs):
	"""organize output such that amino acids are grouped by category"""

	output = []
	for AA in AAs:
		AA_3 = (seq3(AA)).upper()
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
