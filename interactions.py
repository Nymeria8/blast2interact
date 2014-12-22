
from sys import argv

def parse_expression (infile_exp):#recebe ficheiros de expressao
	exp=open(infile_exp)
	dic_exp={}
	for i in exp:
		values=i.split()
		print(values[0])
		dic_exp[values[0]]=values[2]
	exp.close()
	return dic_exp

def read_homology_tab(infile1): #recebe input do blast
	tab=open(infile1)
	homologies={}
	homologies_rev={}
	for i in tab:
		line=i.split("\t")
		if line[0].startswith("#")==False and line[0] not in homologies and float(line[2])>=90.0 and int(line[3])>=100: # identity > 90%; alignment>100aa
			homologies[line[0]]=line[1]
			homologies_rev[line[1]]=line[0]	
	tab.close()
	return homologies, homologies_rev, set(homologies.values()), set(homologies.keys())

homologies, rev, homologies_genes, homologies_vitis=read_homology_tab(argv[1])#recebe o tab do blast


def read_protein_actions(infile2):#recebe string iteraction
	interaction=open(infile2)
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


def generate_edges(graph): #recebe o output do read_protein_actions
	edges = []
	edges_revised=[]
	for node in graph:
		for neighbour in graph[node]:
			if neighbour.endswith("b"):
				edges.append((node, neighbour[:-1],"binding"))
			if neighbour.endswith("r"):
				edges.append((node, neighbour[:-1],"reaction"))
			if neighbour.endswith("c"):
				edges.append((node, neighbour[:-1],"catalysis"))
	for puples in edges:
		empty=(puples[1], puples[0], puples[2])
		if empty not in edges_revised and puples not in edges_revised:
			edges_revised.append(puples)
	print (edges_revised)
	return edges_revised
	
def get_real_names(edges):#recebe o output do generate edges 
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
		

def get_expression(*argv):
	dic_fc={}
	for arg in argv:
		print(arg)
		dic=parse_expression(str(arg).strip("[]").strip("'"))#tem de ir para um dicionario, que depois vai ser adicionado ao outroz
		for value in homologies_vitis:
			if value in dic and value not in dic_fc:
				dic_fc[value]=[dic[value]]
			elif value in dic and value in dic_fc:
				dic_fc[value].append(dic[value])
			elif value not in dic and value in dic_fc:
				dic_fc[value].append("0")
			elif value not in dic and value not in dic_fc:
				dic_fc[value]=["0"]
	return dic_fc

		
def write_expression (outfile2, dic):
	output=open(outfile2, "w")
	for key, value in dic.items():
		output.write(key+"\t")
		for item in value:
			output.write(item+"\t")
		output.write("\n")
	output.close()
    
def write_csv(edges, outfile):
	output=open(outfile, "w")
	empty=()
	for tuples in edges:
		output.write(tuples[0]+"\t"+tuples[1]+"\t"+tuples[2]+"\n")
	output.close()
	



write_csv(get_real_names(generate_edges(read_protein_actions(argv[2]))), argv[3])

write_expression(argv[4], get_expression(argv[5:]))	
