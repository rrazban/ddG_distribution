#!/usr/bin/python

help_msg = 'analyze contact distriubtion of ddGs_res.py residues that reject the null, compare to overall protein contact distribution (Figure 2 in paper)'

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import seaborn as sns

sys.path.append('../')
from Gaussian_mixture import FitGaussianMixture, get_protein_ddG
from ddGs_res import goodness_of_fit, proteins, pretty_proteins

sys.path.append('../utlts')
import residue_info

d_pdb = {'PON': '1v04', 'lipase':'1ex9', '1BTL':'1btl', 'CAII':'1lug', 'DHFR':'1rx2', 'Rnase H':'2rn2', 'Myoglobin':'1a6k', 'Snase':'1stn', 'human Lysozyme':'1rex', 'lysozyme':'1dpx', 'Rnase':'1fs3', 'Barnase':'1a2p', 'AcP':'2acy', 'Ubiqitin':'1ubq', 'Protein G':'2igd', 'Cro_repressor':'1orc'}


def get_contacts(pdb):
	residue = residue_info.Residue('pdbs/{0}.pdb'.format(pdb))
	residues, d_seq = residue.get_residues_sequence()
	contact_matrix, surface_matrix, depth_matrix = residue.get_contact_surface_depth()

	contacts = contact_matrix.sum(axis=1)
	d_contacts = {}
	for res, AA in d_seq.iteritems():
		d_contacts[res] = contacts[res]
	return d_contacts

def make_dataframe(output, rejected_contacts, all_contacts):
	for contact in rejected_contacts:
		output.append([contact, 'reject the null', pretty_proteins[p]])
	for contact in all_contacts:
		output.append([contact, 'all', pretty_proteins[p]])
	return output

def plotout(df_output):
	plt.figure(figsize=(10, 4.5))
	box = sns.boxplot(x="protein", y="contacts", hue="residues", data=df_output)#, palette="PRGn")#, showmeans=True)
	for i,thisbar in enumerate(box.artists):
		if i%2 == 1:
			thisbar.set_hatch('//')
	plt.xlabel('')
	plt.show()


if __name__ == "__main__": 
	data_file = 'Tokuriki_2007.xlsx'
	df = pd.read_excel(data_file)

	output = []
	for p, protein_name in enumerate(proteins):
		protein = df[df['Name']==protein_name]
		ddGs = get_protein_ddG(protein)
		residues_ref = list(protein['No'])

		ks, pvalues = goodness_of_fit(ddGs)
		residues_fail = np.where(np.array(pvalues)<0.05)[0]

		d_contacts = get_contacts(d_pdb[protein_name])
		all_contacts = d_contacts.values()
		rejected_contacts = []	
		for res in residues_fail:
			res_ref = residues_ref[res]
			rejected_contacts.append(d_contacts[res_ref])
		output = make_dataframe(output, rejected_contacts, all_contacts)
#		k, pvalue = stats.ks_2samp(all_contacts, rejected_contacts)	#only CAII significant
	plotout(pd.DataFrame(output, columns = ['contacts', 'residues', 'protein']))
	
