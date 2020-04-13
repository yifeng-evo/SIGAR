#IES file
#bowtie2 insertion, bwa insertion, chimeric supplementary mapping
#bowtie2_insertion_ies=sys.argv[1]
#bwa_insertion_ies=sys.argv[2]
#bwa_supmap_ies=sys.argv[3] #
#cat bowtie2_insertion_ies bwa_insertion_ies bwa_supmap_ies > ies_files
#sort -u ies_files >sorted_ies_files

#using dictionary is better than comparing the adjacent lines


import sys

output_path=sys.argv[1]
ies_output_name=sys.argv[2]
input_file=sys.argv[3]
f=open(output_path+"/process/"+input_file,"r")
f2=open(output_path+"/process/"+input_file+"_simplified_version","w")
f1=open(output_path+ies_output_name,"w")

for line in f:
	f2.write(line.split("\t")[0]+"\t"+line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[7]+"\n")

f2.close()
f.close()

f2=open(output_path+"/process/"+input_file+"_simplified_version","r")

ies_dict={}
for line in f2:
	if line[:-1] not in ies_dict:
		ies_dict[line[:-1]]=1
	else:
		ies_dict[line[:-1]]=ies_dict[line[:-1]]+1

for key in ies_dict:
	f1.write(key+"\t"+str(ies_dict[key])+"\n")

f1.close()
f2.close()









