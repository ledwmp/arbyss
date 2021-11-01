
class residue(object):
	"""Object that holds an amino-acid
	"""
	def __init__(self,atom_list):
		self._atomlist = atom_list
		self._chainID = atom_list[0][4]
		self._resnum = atom_list[0][5]
		self._aa = atom_list[0][3]
		self._atoms = [atom(i) for i in atom_list
		self.oneletter()
	def oneletter(self):
		dict_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q',\
		 			'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T',\
					 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H',\
					  'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A',\
					   'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
		self._aa_one = dict_aa[self._aa]

class pared_residue(residue):
	"""Object that holds a pared-down residue, or beta carbon for non-GLY, alpha carbon for GLY
	"""
	def __init__(self,atom_list):
		super(pared_residue,self).__init__(atom_list)
		if self._aa == "GLY":
			for i in self._atoms:
				if i._atomtype == "CA":
					self._coord = i
		else:
			for i in self._atoms:
				if i._atomtype == "CB":
					self._coord = i

class atom:
	"""Object that holds an atom, including 3D-coordinates
	"""
	def __init__(self,atom_list):
		self._atomID = atom_list[12]
		self._atomtype = atom_list[2]
		self._X = float(atom_list[6])
		self._Y = float(atom_list[7])
		self._Z = float(atom_list[8])
		self._resnum = atom_list[5]
		self._res = atom_list[3]

class chain:
	"""Object that holds a peptide chain
	"""
	def __init__(self,atom_list):
		self._chainID = atom_list[0][4]
		self._atomlist = atom_list[0]
		self._structurelist = atom_list[1]
	def make_residue(self):
		"""Method to take a list of atoms and return a list of residues
		Returns:
			List of residue
		"""
		residue_list = []
		tmp_list = []
		last_res = self._atomlist[0][5]
		for i in self._atomlist:
			if i[5] != last_res:
				residue_list.append(tmp_list)
				tmp_list = []
			else:
				tmp_list.append(i)
			last_res = i[5]
		residue_list.append(tmp_list)
		residue_listfull = [residue(i) for i in residue_list]
		residue_pare = [pared_residue(i) for i in residue_list]
		self._residuelist = residue_listfull
		self._paredlist = residue_pare
		self._chainlen = len(self._residuelist)
	def define_reference(self):
		"""Method to
		Returns:
			List of structure coordinates
		"""
		print(self._structurelist)
		print([i._resnum for i in self._paredlist])
		self._structurecoord = [self._paredlist[i-1]._coord for i in self._structurelist]
	def chunk_out(self,chunk_size):
		"""Method to chunk out sub-peptides
		Args:
			chunk_size: size of peptide surrounding residue
		Returns:
			List of residues with surrounding (n-1)/2
		"""
		tmp_list = []
		self.make_residue()
		self.define_reference()
		print(self._structurecoord)
		len_chain = self._chainlen
		walk = int((chunk_size - 1)/2)
		for i in range(0,len_chain):
			if walk < i < len_chain-walk:
				chunk = self._paredlist[i-walk:i+walk+1]
			elif i < walk:
				chunk = [self._paredlist[0]]*(walk-i) +self._paredlist[0:i+walk+1]
			elif len_chain-walk < i:
				chunk = self._paredlist[i-walk:]+[self._paredlist[-1]]*(len_chain-i+1)
			tmp_list.append(chunk)
		self._chunklist = tmp_list
		return self._chunklist
