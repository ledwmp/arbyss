
class residue(object):
	"""Object that holds an amino-acid
	"""
	def __init__(self,atom_list):
		self.__atomlist = atom_list
		self.__chainID = atom_list[0][4]
		self._aa = atom_list[0][3]
		self._atoms = [atom(i) for i in atom_list]

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
		self.__chainID = atom_list[0][4]
		self.__atomlist = atom_list
	def make_residue(self):
		"""Method to take a list of atoms and return a list of residues
		Returns:
			List of residue
		"""
		residue_list = []
		tmp_list = []
		last_res = self.__atomlist[0][5]
		for i in self.__atomlist:
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
		self._chainlen = len(self.residue_listfull)
	def chunk_out(self,chunk_size):
		len_chain = self._chainlen
		for i in range(0,len_chain):
