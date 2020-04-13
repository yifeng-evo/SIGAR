from Bio import SeqIO
import sys

output_path=sys.argv[1]
maxfreq=int(sys.argv[3])
pointer_file=open(output_path+"/process/2_pointer.txt","w")
f=open(output_path+"/process/sorted_junction_output_script1.txt","r") #sorted_junction file
#####headless, sorted file
line=f.readline()
rname_i_B=line.split("\t")[0]# the name of the first contig
f.close()
f=open(output_path+"/process/sorted_junction_output_script1.txt","r")
lines=f.readlines()
#additional alignment for pointer_flanking MDS
pointer_flanking=int(sys.argv[2])
#threshould_junction_to_end=36, the telo splits shouldn't find a pair. Check the result first
i=0
split_5={}
split_3={}
while i < len(lines):
	entry=lines[i]
	rname_i=entry.split("\t")[0]
	if rname_i==rname_i_B:
		junction_3=entry.split("\t")[3]
		junction_5=entry.split("\t")[2]
		if junction_3!="None":
			qname=entry.split("\t")[1]+"_"+str(i)
			ref_start=int(entry.split("\t")[4])
			ref_end=int(entry.split("\t")[5])
			split_3[qname]=(ref_start,ref_end)
		if junction_5!="None":
			qname=entry.split("\t")[1]+"_"+str(i)
			ref_start=int(entry.split("\t")[4])
			ref_end=int(entry.split("\t")[5])
			split_5[qname]=(ref_start,ref_end)
	else:
		# we should deal with rname_i_B group
		pointer_5_q={}
		pointer_3_q={}
		for qname_3 in split_3:
			for qname_5 in split_5:
				if split_5[qname_5][0]-split_3[qname_3][0]>=pointer_flanking and split_5[qname_5][0]<=split_3[qname_3][1] and split_5[qname_5][1]-split_3[qname_3][1]>=pointer_flanking:				
					pointer_range=(split_5[qname_5][0],split_3[qname_3][1])
					if pointer_range not in pointer_5_q:
						pointer_5_q[pointer_range]=[]
					if qname_5 not in pointer_5_q[pointer_range] and len(pointer_5_q[pointer_range])<=maxfreq:
						pointer_5_q[pointer_range].append(qname_5)
					if pointer_range not in pointer_3_q:
						pointer_3_q[pointer_range]=[]
					if qname_3 not in pointer_3_q[pointer_range] and len(pointer_3_q[pointer_range])<=maxfreq:
						pointer_3_q[pointer_range].append(qname_3)

		pointer_file=open(output_path+"/process/2_pointer.txt","a")
		for key in pointer_5_q:
			pointer_file.write(str(rname_i_B)+"\t"+str(key[0])+"\t"+str(key[1])+"\t"+str(len(pointer_5_q[key]))+"\t"+str(len(pointer_3_q[key]))+"\n")
		split_3={}
		split_5={}
		junction_3=entry.split("\t")[3]
		junction_5=entry.split("\t")[2]
		if junction_3!="None":
			qname=entry.split("\t")[1]+"_"+str(i)
			ref_start=int(entry.split("\t")[4])
			ref_end=int(entry.split("\t")[5])
			split_3[qname]=(ref_start,ref_end)
		if junction_5!="None":
			qname=entry.split("\t")[1]+"_"+str(i)
			ref_start=int(entry.split("\t")[4])
			ref_end=int(entry.split("\t")[5])
			split_5[qname]=(ref_start,ref_end)
	rname_i_B=rname_i
	i=i+1

#######for the last line
if len(split_3)==0 or len(split_5)==0:
	pass
else:
	pointer_5_q={}
	pointer_3_q={}
	for qname_3 in split_3:
		for qname_5 in split_5:
			if split_5[qname_5][0]-split_3[qname_3][0]>=pointer_flanking and split_5[qname_5][0]<=split_3[qname_3][1] and split_5[qname_5][1]-split_3[qname_3][1]>=pointer_flanking:				
				pointer_range=(split_5[qname_5][0],split_3[qname_3][1])
				if pointer_range not in pointer_5_q:
					pointer_5_q[pointer_range]=[]
				if qname_5 not in pointer_5_q[pointer_range] and len(pointer_5_q[pointer_range])<=maxfreq:
					pointer_5_q[pointer_range].append(qname_5)
				if pointer_range not in pointer_3_q:
					pointer_3_q[pointer_range]=[]
				if qname_3 not in pointer_3_q[pointer_range] and len(pointer_3_q[pointer_range])<=maxfreq:
					pointer_3_q[pointer_range].append(qname_3)

	pointer_file=open(output_path+"/process/2_pointer.txt","a")
	for key in pointer_5_q:
		pointer_file.write(str(rname_i_B)+"\t"+str(key[0])+"\t"+str(key[1])+"\t"+str(len(pointer_5_q[key]))+"\t"+str(len(pointer_3_q[key]))+"\n")

f.close()
pointer_file.close()
