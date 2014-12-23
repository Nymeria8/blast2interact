"""
Recives:
-blast output in outfmt-7 against a database of sequences with interactions known,
-a interaction file (gene1 \t gene2 \t interactiontype \n),
-variable number of DESEQ output files of our genes

Returns:
-file with interactions of our sequences, found by homologies:
	-with the names of the sequences in the database
	-with the names of our sequences (default)
-file with genename \t fcvalue1 \t fcvalue2 \t fcvalue3 ....

This files can be imediatly feeded to cytoscape (http://www.cytoscape.org/) for network construction

Usage: python interactions.py blastoutput correspondance_file output_correspondance output_expression DESEQ2_ouput_files[variable]
"""

from sys import argv

def parse_expression (infile_exp):#parse the DESEQ2 output files. output: dictionary name - fold change
	exp=open(infile_exp)
	dic_exp={}
	for i in exp:
		values=i.split()
		print(values[0])
		dic_exp[values[0]]=values[2]
	exp.close()
	return dic_exp

def read_homology_tab(infile1): #parse the tabular file from blast - outfmt=7. returns 1 dictionaries with hit-gene and 2 sets, one with the genes and other with hits
	tab=open(infile1)
	homologies={}
	homologies_rev={}
	for i in tab:
		line=i.split("\t")
		if line[0].startswith("#")==False and line[0] not in homologies and float(line[2])>=90.0 and int(line[3])>=100: # cutoff results: identity > 90%; alignment>100aa
			homologies[line[0]]=line[1]
			homologies_rev[line[1]]=line[0]	
	tab.close()
	return homologies_rev, set(homologies.values()), set(homologies.keys())

rev, homologies_genes, homologies_vitis=read_homology_tab(argv[1])#returns the read_homology_tab, as global, so it can be used by several other fuctions


def read_protein_actions(infile2):#parses the interaction file. http://string-db.org/ file or any other file with element one {tab} element two {tab} type of interaction
	interaction=open(infile2)#returns a dictionary with the interactions of the homologies found
	pairs={}
	for i in interaction:
		line=i.split("\t")
		if line[0].startswith("item")==False and line[0] not in pairs and line[0] in homologies_genes and line[1] in homologies_genes:
			if line[2]=="binding":
				pairs[line[0]]=set([line[1]+"b"])
			if line[2]=="reaction":
				pairs[line[0]]=set([line[1]+"r"])
			if line[2]=="catalysis":
				pairs[line[0]]=set([line[1]+"c"])
		elif line[0].startswith("item")==False and line[0] in pairs and line[1] in homologies_genes:
			if line[2]=="binding":
				pairs[line[0]].add(line[1]+"b")
			if line[2]=="reaction":
				pairs[line[0]].add(line[1]+"r")
			if line[2]=="catalysis":
				pairs[line[0]].add(line[1]+"c")
	interaction.close()
	return pairs


def generate_edges(graph): #parses the output of read_protein_actions, and identifies the homologies of our sequences. returns a list of tuples with the interactions of the sequences by the homologies, not the genes
	edges = []
	edges_revised=[]
	for node in graph:
		for neighbour in graph[node]:
			if neighbour.endswith("b"):
				edges.append((node, neighbour[:-1],"binding"))#type of interaction
			if neighbour.endswith("r"):
				edges.append((node, neighbour[:-1],"reaction"))#type of interaction
			if neighbour.endswith("c"):
				edges.append((node, neighbour[:-1],"catalysis"))#type of interaction
	for puples in edges:
		empty=(puples[1], puples[0], puples[2])
		if empty not in edges_revised and puples not in edges_revised:
			edges_revised.append(puples)
	print (edges_revised)#because of the different types of interaction, a pair of genes can apear several times but with a different interaction
	return edges_revised
	
def get_real_names(edges):#uses the output of generate_edges. retrives the  same thing of generate_edges, but with the identifiers of our genes
	l=[]
	final_list=[]
	for tuples in edges:
		for element in tuples:
			if element.startswith("binding") or element.startswith("reac") or element.startswith("cat"):
				l.append(element)
				tup=tuple(l)
				final_list.append(tup)
				l=[]
			else:
				if element in rev:
					l.append(rev[element])
	return final_list
		

def get_expression(*argv):#call parse_expression to read the DESEQ output and select the ones correspondant to the used genes. if the gene dont apear in the output of DESEQ, a value of 0 is automaticly gived
	dic_fc={}#recives as many as needed DESEQ outputs
	for arg in argv:
		print(arg)
		dic=parse_expression(str(arg).strip("[]").strip("'"))#tem de ir para um dicionario, que depois vai ser adicionado ao outroz
		for value in homologies_vitis:
			if value in dic and value not in dic_fc:
				dic_fc[value]=[dic[value]]
			elif value in dic and value in dic_fc:
				dic_fc[value].append(dic[value])
			elif value not in dic and value in dic_fc:
				dic_fc[value].append("0")#absent from the DESEQ output - value 0
			elif value not in dic and value not in dic_fc:
				dic_fc[value]=["0"]#absent from the DESEQ output - value 0
	return dic_fc#retrives a dictionary: gene:[value1, value2, value3...]

		
def write_expression (outfile2, dic):#recives the dictionary from get expression, and write a file with gene_name {tab} value1 {tab} value2 ....
	output=open(outfile2, "w")
	for key, value in dic.items():
		output.write(key+"\t")
		for item in value:
			output.write(item+"\t")
		output.write("\n")
	output.close()
    
def write_csv(edges, outfile):#recives the get_real_names output (default) or read_protein_actions and write a file with name1 {tab} name2 {tab} interaction {enter}
	output=open(outfile, "w")
	empty=()
	for tuples in edges:
		output.write(tuples[0]+"\t"+tuples[1]+"\t"+tuples[2]+"\n")
	output.close()
	



write_csv(get_real_names(generate_edges(read_protein_actions(argv[2]))), argv[3])# change here to change from the bd names or our sequence names

write_expression(argv[4], get_expression(argv[5:]))	
