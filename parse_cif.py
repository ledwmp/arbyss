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

def whole_chain(tmp_list):
	chain_list = []
	if tmp_list[0] == "loop_\n":
		keys = [i.strip() for i in tmp_list if i[0] == "_"]
		values = [[j for j in i.split(" ") if j != ""] for i in tmp_list[1:] if i[0] != "_"]
	labels = ["_entity_poly_seq.entity_id","_entity_poly_seq.num",\
			"_entity_poly_seq.mon_id","_entity_poly_seq.hetero"
			]
	labels = [keys.index(i) for i in labels]
	aa_index = keys.index("_entity_poly_seq.mon_id")
	for i in values:
		if i[0] == "1":
			chain_list.append(i[aa_index])
	return chain_list


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

def parse_helices(tmp_list):
	structure_list = defaultdict(list)
	if tmp_list[0] == "loop_\n":
		keys = [i.strip() for i in tmp_list if i[0] == "_"]
		values = [[j for j in i.split(" ") if j != ""] for i in tmp_list[1:] if i[0] != "_"]
	labels = ["_struct_conf.conf_type_id","_struct_conf.id","_struct_conf.pdbx_PDB_helix_id",\
		"_struct_conf.beg_label_comp_id","_struct_conf.beg_label_asym_id","_struct_conf.beg_label_seq_id",\
		"_struct_conf.pdbx_beg_PDB_ins_code","_struct_conf.end_label_comp_id","_struct_conf.end_label_asym_id",\
		"_struct_conf.end_label_seq_id","_struct_conf.pdbx_end_PDB_ins_code","_struct_conf.beg_auth_comp_id",\
		"_struct_conf.beg_auth_asym_id","_struct_conf.beg_auth_seq_id","_struct_conf.end_auth_comp_id",\
		"_struct_conf.end_auth_asym_id","_struct_conf.end_auth_seq_id","_struct_conf.pdbx_PDB_helix_class",\
		"_struct_conf.details","_struct_conf.pdbx_PDB_helix_length"
		]
	labels = [keys.index(i) for i in labels]
	chain_index = keys.index('_struct_conf.beg_label_asym_id')
	aa_start_index = keys.index("_struct_conf.beg_label_seq_id")
	aa_end_index = keys.index("_struct_conf.end_label_seq_id")
	for i in values:
		if i[0] == "HELX_P":
			structure_list[i[chain_index]].append(i[aa_start_index])
			structure_list[i[chain_index]].append(i[aa_end_index])
	return structure_list

def parse_sheets(tmp_list):
	structure_list = defaultdict(list)
	if tmp_list[0] == "loop_\n":
		keys = [i.strip() for i in tmp_list if i[0] == "_"]
		values = [[j for j in i.split(" ") if j != ""] for i in tmp_list[1:] if i[0] != "_"]
	labels = ["_struct_sheet_range.sheet_id","_struct_sheet_range.id",\
		"_struct_sheet_range.beg_label_comp_id","_struct_sheet_range.beg_label_asym_id",\
		"_struct_sheet_range.beg_label_seq_id","_struct_sheet_range.pdbx_beg_PDB_ins_code",\
		"_struct_sheet_range.end_label_comp_id","_struct_sheet_range.end_label_asym_id",\
		"_struct_sheet_range.end_label_seq_id","_struct_sheet_range.pdbx_end_PDB_ins_code",\
		"_struct_sheet_range.beg_auth_comp_id","_struct_sheet_range.beg_auth_asym_id",\
		"_struct_sheet_range.beg_auth_seq_id","_struct_sheet_range.end_auth_comp_id",\
		"_struct_sheet_range.end_auth_asym_id","_struct_sheet_range.end_auth_seq_id"
		]
	labels = [keys.index(i) for i in labels]
	chain_index = keys.index('_struct_sheet_range.beg_label_asym_id')
	aa_start_index = keys.index("_struct_sheet_range.beg_label_seq_id")
	aa_end_index = keys.index("_struct_sheet_range.end_label_seq_id")
	for i in values:
		if i[0] == "A":
			structure_list[i[chain_index]].append(i[aa_start_index])
			structure_list[i[chain_index]].append(i[aa_end_index])
	return structure_list

def list_of_atoms(in_file,chain_id):
	pdb_list = read_pdb(in_file)
	whole_chain([i for i in pdb_list if "_entity_poly_seq.entity_id \n" in i][0])
	list_of_atoms = parse_atoms([i for i in pdb_list if "_atom_site.group_PDB \n" in i][0])
	list_of_helices = parse_helices([i for i in pdb_list if "_struct_conf.conf_type_id \n" in i][0])
	list_of_sheets = parse_sheets([i for i in pdb_list if "_struct_sheet_range.sheet_id \n" in i][0])
	sheets_helices = defaultdict(list)
	for k,v in list_of_helices.items():
		for i in v:
			sheets_helices[k].append(i)
	for k,v in list_of_sheets.items():
		for i in v:
			sheets_helices[k].append(i)
	return list_of_atoms[chain_id],sheets_helices[chain_id]
