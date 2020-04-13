import sys
output_path=sys.argv[1]
output_name=sys.argv[2]

f3=open(output_path+"/process/sorted_junction_output_script1.txt","r")#"0801_junction_BWA_MIC_to_MAC_Ewoo_bigthan500_telo1or2.txt"
output=open(output_path+"/results/junction_summary_"+output_name,"w")

jct={}

for line in f3:
	junction_5=line.split("\t")[2]
	junction_3=line.split("\t")[3]
	rname=line.split("\t")[0]
	if rname not in jct:
		jct[rname]=[junction_3,junction_5]
	else:
		jct[rname].append(junction_3)
		jct[rname].append(junction_5)

for key in jct:
	l=jct[key]
	jct[key]=dict((x,l.count(x)) for x in set(l))

for rname in jct:
	for junction in jct[rname]:
		if junction!="None":
			output.write(rname+"\t"+junction+"\t"+str(jct[rname][junction])+"\n")

output.close()
f3.close()



