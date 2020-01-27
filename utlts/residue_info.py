import os, sys
from Bio.PDB import PDBParser, NeighborSearch 
from Bio import PDB
import numpy as np

import pdb_info


class Residue(pdb_info.Protein):
	""" 
	This is a class that calculates contact matrices
	of protein residues 
      
	Attributes: 
		contact_matrix (array-like): matrix of 0 and 1s 
			denoting whether two residues are in contact 
	"""

	def __init__(self, pdb_path):
		super(Residue, self).__init__(pdb_path)
		self.max_residue = max(self.residues)+1
		self.pdb_path = pdb_path

		self.atom_list = []
		self.atoms()
		self.contact_matrix = np.zeros((self.max_residue, self.max_residue), dtype=int)

	def atoms(self):
		"""obtain all atoms in the protein"""

		for residue in self.structure.get_residues():
			if PDB.is_aa(residue, standard=True):
				self.atom_list.extend(residue.get_atoms())

	def contact_pairs(self):
		"""obtain all residue pairs with atoms in close proximity"""

		radius = 4.5
		ns = NeighborSearch(self.atom_list)
		self.contact_pairs = ns.search_all(radius, level = 'R')	#residues returned

	def contact(self):
		"""create contact matrix. omit adjacent residues"""

		self.contact_pairs()
		M = np.zeros((self.max_residue, self.max_residue), dtype=int)

		for contact in self.contact_pairs:
			res1 = contact[0].id[1]
			res2 = contact[1].id[1]
			if not abs(res1 - res2) in [1, 0]:	#no nearest neighbors
				if self.contact_matrix[res1][res2]==1:
					print 'already'
				M[res1][res2] = 1
		self.contact_matrix = M + M.T
	
	def get_contact_matrix(self):
		self.contact()
		return self.contact_matrix
