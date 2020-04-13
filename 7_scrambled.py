#scrambled\
import sys
#add qname in the scrmabled _by_deletion part

#then check if the mds junction are already in the non-scrambled pointer junctions, if not, keep them
output_path=sys.argv[1]# the directory to the output file
output_name=sys.argv[2]

chimeric_scrambled_file=open(output_path+"/process/scrambling_chimeric.txt","r")#"scrambling_nosecondary_q10_BWA_MIC_to_MAC_Ewoo.txt"
nonscrambled_revised_format=open(output_path+"/process/mds_ies_bwa_same_contig.txt","r")#MDS_IES_nosecondary_q10_BWA_MIC_to_MAC_Ewoo_sam_same_contig_revised_format.txt
mds_3_or_more=open(output_path+"/process/3_mds_blocks.txt","r")#"nosecondary_q10_BWA_MIC_to_MAC_Ewoo_sam_same_contig_3_mds_or_more.txt"
output=open(output_path+"/results/scrambled_summary_"+output_name,"w")

mds_3_or_more_dict={}
for line in mds_3_or_more:
	rname=line.split("\t")[0]
	qname=line.split("\t")[1]
	if rname not in mds_3_or_more_dict:
		mds_3_or_more_dict[rname]=[qname]
	else:
		mds_3_or_more_dict[rname].append(qname)

non_scrambled_dict={}
for line in nonscrambled_revised_format:
	rname=line.split("\t")[0]
	qname=line.split("\t")[1]
	if rname in mds_3_or_more_dict and qname in mds_3_or_more_dict[rname]:
		pointer_start=int(line.split("\t")[2])
		pointer_end=int(line.split("\t")[3])
		if rname not in non_scrambled_dict:
			non_scrambled_dict[rname]={qname:[(pointer_start,pointer_end)]}
		else:
			if qname not in non_scrambled_dict[rname]:
				non_scrambled_dict[rname][qname]=[(pointer_start,pointer_end)]
			else:
				non_scrambled_dict[rname][qname].append((pointer_start,pointer_end))

for line in chimeric_scrambled_file:
	rname=line.split("\t")[0]
	qname=line.split("\t")[5]
	mds_1_end_boundary=0
	mds_2_start_boundary=0
	if rname in non_scrambled_dict:
		if qname in non_scrambled_dict[rname]:
			mds_1_end=int(line.split("\t")[2])
			mds_2_start=int(line.split("\t")[3])
			for entry in non_scrambled_dict[rname][qname]:
				pointer_start=entry[0]
				pointer_end=entry[1]
				if mds_1_end==pointer_end:
					mds_1_end_boundary=1
				if mds_2_start==pointer_start:
					mds_2_start_boundary=1
	if mds_2_start_boundary==0 and mds_1_end_boundary==0:
		output.write(line)

chimeric_scrambled_file.close()
nonscrambled_revised_format.close()
mds_3_or_more.close()
output.close()







