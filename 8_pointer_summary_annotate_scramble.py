import sys

output_path=sys.argv[1]# the directory to the output file
output_name=sys.argv[2]

f1=open(output_path+"/results/scrambled_summary_"+output_name,"r")#"scrambled_summary_Ewoo_092118_without_deletion_inferred"
f2=open(output_path+"/process/pointer_summary_"+output_name,"r")#"pointer_summary_Ewoo_split_read_mapping_20180918.txt"
output2=open(output_path+"/results/pointer_summary_"+output_name+"_with_scrambled.txt","w")

scrambled_jct={}
for line in f1:
	rname=line.split("\t")[0]
	mds1_end=int(line.split("\t")[2]) # these junctions are supposed to be junctions from split reads
	mds2_start=int(line.split("\t")[3])
	if rname not in scrambled_jct:
		scrambled_jct[rname]=[mds1_end,mds2_start]
	else:
		scrambled_jct[rname].append(mds1_end)
		scrambled_jct[rname].append(mds2_start)

for line in f2:
	rname=line.split("\t")[0]
	if rname in scrambled_jct:
		jct_1=int(line.split("\t")[1])
		jct_2=int(line.split("\t")[2])
		if jct_1 in scrambled_jct[rname] or jct_2 in scrambled_jct[rname]:
			if "IES" not in line: # if there's non scrambled detected, don't report non-scrambled
				output2.write(line[:-1]+"\t"+"scrambled\n")
			else:
				output2.write(line)
		else:
			output2.write(line)
	else:
		output2.write(line)

output2.close()
f1.close()
f2.close()

	





