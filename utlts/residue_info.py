help_msg = 'inherit contact matrix of a PDB'

import os, sys
from Bio.PDB import PDBParser, NeighborSearch 
from Bio import PDB
import numpy as np

import pdb_info


class Residue(pdb_info.Protein):
	def __init__(self, pdb_path):
		super(Residue, self).__init__(pdb_path)
		self.max_residue = max(self.residues)+1
		self.pdb_path = pdb_path

		self.atom_list = []
		self.atoms()
		self.contact_matrix = np.zeros((self.max_residue, self.max_residue), dtype=int)
		self.depth_matrix = np.zeros(self.max_residue)

	def atoms(self):
		for residue in self.structure.get_residues():
			if PDB.is_aa(residue, standard=True):   #only consider actual residues
				self.atom_list.extend(residue.get_atoms())

	def contact_pairs(self):
		radius = 4.5
		ns = NeighborSearch(self.atom_list)
		self.contact_pairs = ns.search_all(radius, level = 'R')	#residues returned

	def contact(self):
		self.contact_pairs()
		M = np.zeros((self.max_residue, self.max_residue), dtype=int)	#in case radius changed

		for contact in self.contact_pairs:
			res1 = contact[0].id[1]
			res2 = contact[1].id[1]
			if not abs(res1 - res2) in [1, 0]:	#no nearest neighbors
				if self.contact_matrix[res1][res2]==1:
					print 'already'
				M[res1][res2] = 1
		self.contact_matrix = M + M.T	#cant set equal to self or else messes up
	
	def get_contact_matrix(self):
		self.contact()
		return self.contact_matrix
