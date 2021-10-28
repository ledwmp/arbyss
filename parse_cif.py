from collections import defaultdict

def read_pdb(pdb_file):
	tmp_list = []
	with open(pdb_file) as r:
		split_list = []
		for line in r:
			if line[0] == "#":
				tmp_list.append(split_list)
				split_list = []
			else:
				split_list.append(line)
	r.close()
	return tmp_list

def parse_atoms(tmp_list):
	chain_list = defaultdict(list)
	#convert_dict = defaultdict(list)
	if tmp_list[0] == "loop_\n":
		keys = [i.strip() for i in tmp_list if i[0] == "_"]
		values = [[j for j in i.split(" ") if j != ""] for i in tmp_list[1:] if i[0] != "_"]
	labels = ['_atom_site.group_PDB','_atom_site.id','_atom_site.label_atom_id',\
		'_atom_site.label_comp_id','_atom_site.auth_asym_id','_atom_site.auth_seq_id', \
		'_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z',\
		'_atom_site.occupancy','_atom_site.B_iso_or_equiv','_atom_site.label_asym_id',\
		'_atom_site.type_symbol'
		]
	labels = [keys.index(i) for i in labels]
	chain_index = keys.index('_atom_site.auth_asym_id')
	ref_index = keys.index("_atom_site.label_seq_id")
	model_index = keys.index("_atom_site.auth_seq_id")
	for i in values:
	 	if i[0] == "ATOM":
			chain_list[i[chain_index]].append([i[j] for j in labels])
			#if i[model_index] != "." and i[ref_index] != ".":
				#convert_dict[i[chain_index]].append((i[model_index],i[ref_index]))
	return chain_list#,convert_dict

def list_of_atoms(in_file,chain_id):
	pdb_list = read_pdb(in_file)
	list_of_atoms = parse_atoms([i for i in pdb_list if "_atom_site.group_PDB \n" in i][0])
	return list_of_atoms[chain_id]

