
class residue(object):
	"""Object that holds an amino-acid
	"""
	def __init__(self,rawnum,atom_list):
		self._atomlist = atom_list
		self._chainID = atom_list[0][4]
		self._resnum = atom_list[0][5]
		self._resraw = rawnum
		self._aa = atom_list[0][3]
		self._atoms = [atom(i) for i in atom_list]
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
	def __init__(self,rawnum,atom_list):
		super(pared_residue,self).__init__(rawnum,atom_list)
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

class alignment_block:
	"""Object that holds a substring to be used for structural alignment
	"""
	def __init__(self,res_list,aa_one,resraw,walk_outer):
		self._res_list = res_list #list of pared residues in substring
		self._aa_one = aa_one #aa
		self._resraw = resraw #raw position in peptide
		self._walk_outer = walk_outer #size of substring for alignments


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
		residue_listfull = [residue(*i) for i in enumerate(residue_list)]
		residue_pare = [pared_residue(*i) for i in enumerate(residue_list)]
		self._residuelist = residue_listfull
		self._paredlist = residue_pare
		self._chainlen = len(self._residuelist)

	def chunk_out(self,chunk_size_outer):
		"""Method to chunk out sub-peptides
		Args:
			chunk_size: size of peptide surrounding residue
		Returns:
			List of residues with surrounding (n-1)/2
		"""
		tmp_list = []
		self.make_residue()
		len_chain = self._chainlen
		walk_outer = int((chunk_size_outer-1)/2)
		self._walk_outer = walk_outer
		for i in range(0,len_chain):
			if walk_outer <= i < len_chain-walk_outer:
				chunk = self._paredlist[i-walk_outer:i+walk_outer+1]
			elif i < walk_outer:
				chunk = self._paredlist[0:i+walk_outer+1]
			elif len_chain-walk_outer <= i:
				chunk = self._paredlist[i-walk_outer:]
			tmp_list.append(alignment_block(\
							chunk,\
							self._paredlist[i]._aa_one,\
							self._paredlist[i]._resraw,\
							walk_outer,\
							))
		self._chunklist = tmp_list
		self._chunklist_dict = {i._resraw:i for i in self._chunklist}
		return self._chunklist
