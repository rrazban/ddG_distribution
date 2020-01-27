import os
from Bio.PDB import PDBParser, Residue, Polypeptide
from Bio import PDB


class Protein(object):
	""" 
	This is a class that obtains protein sequence 
	and structure information 
      
	Attributes: 
		d_seq (dict of {int: str): connects PDB 
			residue number with its aminoacid type
	"""

	def __init__(self, pdb_path):
		self.structure = PDBParser().get_structure("", pdb_path)
		self.residues = []
		self.d_sequence = {}

		self.parse_structure()

	def parse_structure(self):
		for residue in self.structure.get_residues():
			if PDB.is_aa(residue, standard=True):	#only consider standard 20 residues
				res = residue.id[1]
				if res not in self.residues:	#dont doublecount mutated residues	(ex. 1ORC)	
					self.residues.append(res)
					self.d_sequence[res] = Polypeptide.three_to_one(Residue.Residue.get_resname(residue))

	def get_residues_sequence(self):
		return self.residues, self.d_sequence
