#!/usr/bin/python3.5
 
import numpy as np
import re
import os
import enum
import sys
from Bio import SeqIO

arg = sys.argv
print(arg)
inputFile = arg[1]
originFile = arg[2]
outputFile1 = arg[3]
outputFile2 = arg[4]


#------Extract gene name, contig name, CIGAR put it on a tab------#
def get_data(file):
	tab = []
	with open(file, "r") as f1 :
		for line in f1:
			if re.search('^NODE', line):
				tabline = line.split("\t")
				if re.search('^ARO', tabline[2]):
					aro = tabline[2].split("|")
					tab.append(tabline[0]+"\t"+aro[2]+"\t"+tabline[5])
	f1.close()
	return tab 	


def get_arg_length (file): 
	dic_arg = {}
	for i in SeqIO.parse(file, 'fasta'):
		name = i.id
		seq = i.seq
		seqLen = len(i)
		aro=name.split("|")
		dic_arg[aro[2]] = seqLen
	return dic_arg	

# Enum for size units
class SIZE_UNIT(enum.Enum):
   BYTES = 1
   KB = 2
   MB = 3
   GB = 4
   
#------convert file size------#
def convert_unit(size_in_bytes, unit):
   if unit == SIZE_UNIT.KB:
       return size_in_bytes/1024
   elif unit == SIZE_UNIT.MB:
       return size_in_bytes/(1024*1024)
   elif unit == SIZE_UNIT.GB:
       return size_in_bytes/(1024*1024*1024)
   else:
       return size_in_bytes
       	
def get_file_size(file_name, size_unit=SIZE_UNIT.BYTES):
   size = os.path.getsize(file_name)
   return convert_unit(size, size_unit)

def get_arg_length_fasta (file): 
	dic_arg = {}
	for i in SeqIO.parse(file, 'fasta'):
		name = i.id
		seq = i.seq
		seqLen = len(i)
		dic_arg[name] = seqLen
	return dic_arg	

def get_abundance(sample_size, read_size ,ARGsize, nb_mapped_read):
	ARGsize = int(ARGsize)
	nb_mapped_read = int(nb_mapped_read)
	read_size = int(read_size)
	#print(sample_size, read_size ,ARGsize, nb_mapped_read)
	abund = (nb_mapped_read*read_size/ARGsize)/sample_size
	return abund
	
	
def main():
#------STEP1------#
	#recup NODE, ARG,  CIGAR, taille NODE
	tab = []
	tab = get_data(inputFile)
	
	#recup les tailles de la sequence des ARG
	dic_arg = {}
	dic_arg = get_arg_length("../bdd/wildcard_database_v3.0.1.fasta")

	#recup taille du fichier
	file_path = originFile
	sample_size = get_file_size(file_path, SIZE_UNIT.GB)
	print("size of "+file_path+"\t"+":"+"\t"+str(sample_size)+"\t"+"Gb.")
	
	dic_node = {}
	dic_node = get_arg_length_fasta(originFile)

	#merge toutes les data : ARG, NODE, CIGAR, taille NODE, taille ARG
	tab_final = []
	list_arg = []
	for i in tab : 
		data = i.split("\t")
		if data[1] in dic_arg.keys():
			if data[0] in dic_node.keys():
				list_arg.append(data[1])
				tab_final.append(data[1]+"\t"+data[0]+"\t"+data[2]+"\t"+str(dic_node[data[0]])+"\t"+str(dic_arg[data[1]]))
	#print (tab_final)
	
	#recup nb reads qui mappe sur ARG
	di = dict((x,list_arg.count(x)) for x in set(list_arg))
	#print (di)
	
	#initialisation des abundance par ARG a 0
	dic_abundByARG = {}
	for arg in di :
		dic_abundByARG[arg] = 0

	tabFINAL = []
	for i in tab_final :
		data = i.split("\t")
		ARGname = data[0]
		NODEname = data[1]
		CIGAR = data[2]
		NODEsize = data[3]
		ARGsize =data[4]
		nbMappedContig = di[data[0]]
		abund = get_abundance(sample_size, NODEsize, ARGsize, nbMappedContig)
		if ARGname in dic_abundByARG :
			dic_abundByARG[ARGname] = float(dic_abundByARG[ARGname]) + float(abund)
		tabFINAL.append(ARGname+"\t"+NODEname+"\t"+CIGAR+"\t"+str(NODEsize)+"\t"+str(ARGsize)+"\t"+str(nbMappedContig)+"\t"+str(abund))
	#print(tabFINAL)

	with open(outputFile1, "w") as fout:  
		fout.write("ARGname"+"\t"+"NODEname"+"\t"+"CIGAR"+"\t"+"NODEsize"+"\t"+"ARGsize"+"\t"+"nbMappedContig"+"\t"+"ABUNDANCE"+"\t"+"ABUNDANCE TOTALE"+"\n")
		for i in tabFINAL: 
			data = i.split("\t")
			if data[0] in dic_abundByARG :
				fout.write(str(i)+"\t"+str(dic_abundByARG[data[0]])+"\n")
	fout.close()

	
	with open(outputFile2, "w") as fabund:  
		fabund.write("gene"+"\t"+"nbmappedContig"+"\t"+"abundance"+"\n")
		for arg in di :
			if arg in dic_abundByARG:
				fabund.write(str(arg)+"\t"+str(di[arg])+"\t"+str(dic_abundByARG[arg])+"\n")
	fabund.close()
	
	
if __name__ == '__main__':
	main()