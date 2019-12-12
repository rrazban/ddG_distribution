help_msg = 'inherit residue, sequence, chain info about PDB'	

import os
from Bio.PDB import PDBParser, Residue, Polypeptide
from Bio import PDB


AAs = "LFIMVWCYHATGPRQSNEDK"	#not used here but other scripts that import this module

class Protein(object):
	def __init__(self, pdb_path):
#		self.pdb = os.path.basename(pdb_path)[:-4]
#		self.chain = self.pdb[self.pdb.index('-')+1:self.pdb.index('_')]	#specific for Repair jobs

		self.structure = PDBParser().get_structure("", pdb_path)
		self.residues = []		#only includes those residues captured by Xray crystallography
		self.d_sequence = {}

		self.parse_structure()

	def parse_structure(self):
		for residue in self.structure.get_residues():
			if PDB.is_aa(residue, standard=True):	#only consider actual residues	#only consider standard 20
				res = residue.id[1]
				if res not in self.residues:	#dont doublecount mutated residues	(ex. 1ORC)	
					self.residues.append(res)
					self.d_sequence[res] = Polypeptide.three_to_one(Residue.Residue.get_resname(residue))

	def get_residues_sequence(self):
		return self.residues, self.d_sequence#, self.chain