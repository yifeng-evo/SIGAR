#cat ies by insertions
#sort -u ies_by_insertions 
# IES pipe for ies_by_insertion and ies_by _chimer_reads
#report number of reads supporting the pointer; and if there are ies detected



import sys
from Bio import SeqIO
output_path=sys.argv[1]
ies_summary=sys.argv[2]
output=sys.argv[3]
MAC_genome=sys.argv[4]

f1=open(output_path+"/process/2_pointer.txt","r")
f3=open(output_path+"/results/"+ies_summary,"r")
ies_summary_list=[]


for line in f3:
	ies_summary_list.append(line.split("\t")[0:3])

output_file_withies=open(output_path+'/process/'+output,"w")
genome_dict={}
for record in SeqIO.parse(MAC_genome,"fasta"):
	genome_dict[record.id]=record.seq

for line in f1:
	rname=line.split("\t")[0]
	pointer_start=int(line.split("\t")[1])
	pointer_end=int(line.split("\t")[2])
	if line.split("\t")[0:3] in ies_summary_list:
		output_file_withies.write(line[:-1]+"\t"+str(genome_dict[rname][pointer_start:pointer_end])+"\tIES available\n")
	else:
		output_file_withies.write(line[:-1]+"\t"+str(genome_dict[rname][pointer_start:pointer_end])+"\n")
f1.close()
output_file_withies.close()

