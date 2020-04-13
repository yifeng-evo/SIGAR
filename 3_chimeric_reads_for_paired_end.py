import re
import math
import sys
#sys.argv[1]:outputpath

output_path=sys.argv[1]
sam=open(output_path+sys.argv[2],"r")#filtered bwa sam


mapping={}

def end_on_contig(CIGAR, start, strand,mapq):
	matchM = re.findall(r'(\d+)M', CIGAR)
	intlistM = [int(x) for x in matchM]
	sumM=sum(intlistM)
	matchD = re.findall(r'(\d+)D', CIGAR)
	intlistD = [int(x) for x in matchD]
	sumD=sum(intlistD)
	end=sumM+sumD+start-1


	matchI = re.findall(r'(\d+)I', CIGAR)
	intlistI = [int(x) for x in matchI]
	sumI=sum(intlistI)
	soft_clip_5=0
	soft_clip_3=0
	hard_clip_5=0
	hard_clip_3=0
	if "S" in CIGAR:
		if CIGAR.count("S")==1:
			if CIGAR.split("S")[0].isdigit()==True: #soft clip at 5' end
				soft_clip_5=int(CIGAR.split("S")[0])
			else:
				soft_clip_3=int(re.findall(r'(\d+)S', CIGAR)[0])
		elif CIGAR.count("S")==2:
			soft_clip_5=int(re.findall(r'(\d+)S', CIGAR)[0])
			soft_clip_3=int(re.findall(r'(\d+)S', CIGAR)[1])

	if "H" in CIGAR:
		if CIGAR.count("H")==1:
			if CIGAR.split("H")[0].isdigit()==True: #soft clip at 5' end
				hard_clip_5=int(CIGAR.split("H")[0])
			else:
				hard_clip_3=int(re.findall(r'(\d+)H', CIGAR)[0])
		elif CIGAR.count("H")==2:
			hard_clip_5=int(re.findall(r'(\d+)H', CIGAR)[0])
			hard_clip_3=int(re.findall(r'(\d+)H', CIGAR)[1])

	clip_3=soft_clip_3+hard_clip_3
	clip_5=hard_clip_5+soft_clip_5

	if strand==1: #plus strand
		query_start=clip_5
	else:
		query_start=clip_3
	query_end=query_start+sumM+sumI


	return [(start-1,end),(query_start,query_end),mapq,strand]


def strand (flag):
	plus=1 #plus=1: true; plus=0: reverse. defalut: plus strand
	a=10
	while a >4 and flag >0:
		a=int(math.log(flag,2))
		if a == 4:
			plus=0
		flag=flag-2**a

	return plus

def pair(flag):
	pair=0 #first=1: first in pair; first=2: second in pair
	a=10
	while a >6 and flag >0:
		a=int(math.log(flag,2))
		if a == 6:
			pair=1
		if a==7:
			pair=2
		flag=flag-2**a

	return pair

def uniq(input):
	output = []
	for x in input:
		if x not in output:
			output.append(x)
	return output




for line in sam:
	if "SA:Z:" in line:
		rname=line.split("\t")[2]
		SAsection=line.split("\t")[15]
		if rname in SAsection:
			qname=line.split("\t")[0]
			flag=int(line.split("\t")[1])
			CIGAR=line.split("\t")[5]
			mapq=int(line.split("\t")[4])
			start=int(line.split("\t")[3])
			#if pair(flag)!=0:
				#qname=qname+"/"+str(pair(flag))
			if rname in mapping:
				if qname not in mapping[rname]:
					mapping[rname][qname]=[]
			else:
				mapping[rname]={}
				mapping[rname][qname]=[] 

			span=end_on_contig(CIGAR,start,strand(flag),mapq) # don't consider mapq>=20 because this is already filtered
			if span not in mapping[rname][qname]:
				mapping[rname][qname].append(span)
		#mapping is to record the supplementary records

#chimeric_reads.write("rname\tqname\tr_span\tq_span\tmapq\n")
for rname in mapping:
	for qname in mapping[rname]:
		mapping[rname][qname]=uniq(mapping[rname][qname])
#		for span in mapping[rname][qname]:
#			if span[2]>=10:
#				chimeric_reads.write(rname+"\t"+qname+"\t"+str(span[0])+"\t"+str(span[1])+"\t"+str(span[2])+"\t"+str(span[3])+"\n")



sam.close()
#chimeric_reads.close()


MDS_IES=open(output_path+"/process/mds_ies_bwa_same_contig.txt","w")#MDS_IES_nosecondary_q10_BWA_MIC_to_MAC_Ewoo_sam_same_contig_revised_format.txt
scrambling=open(output_path+"/process/scrambling_chimeric.txt","w")#"scrambling_nosecondary_q10_BWA_MIC_to_MAC_Ewoo.txt"
telocontig=open(sys.argv[3],"r")#"contigs_name_bigthan500_telo1or2_Ewoo_MAC_v4.3.txt"
MDS_3_blocks=open(output_path+"/process/3_mds_blocks.txt","w")#"nosecondary_q10_BWA_MIC_to_MAC_Ewoo_sam_same_contig_3_mds_or_more.txt"
telomeric_contig=[]
for line in telocontig:
	telomeric_contig.append(line.split("\t")[0])


for rname in mapping:
	if rname in telomeric_contig:
		for qname in mapping[rname]:
			number_records=len(mapping[rname][qname])
			if number_records >=3:
				n=0
				while n < number_records:
					MDS_3_blocks.write(rname+"\t"+qname+"\t"+str(mapping[rname][qname][n][0])+"\t"+str(mapping[rname][qname][n][1])+"\t"+str(mapping[rname][qname][n][2])+"\t"+str(mapping[rname][qname][n][3])+"\n")
					n=n+1
			i=0
			while i < number_records:
				k=i+1
				while k< number_records:
					if mapping[rname][qname][k][3]==mapping[rname][qname][i][3] and mapping[rname][qname][k][3]==1: # same plus strand
						if mapping[rname][qname][k][0][0]< mapping[rname][qname][i][0][0]<=mapping[rname][qname][k][0][1] and mapping[rname][qname][k][0][1]<mapping[rname][qname][i][0][1] and mapping[rname][qname][k][1][1]<mapping[rname][qname][i][1][0]:#ith record rstart > i+1th record rstart
							#rname/pointer_Start/pointer_end/qname/IES_start/IESend
							MDS_IES.write(rname+"\t"+qname+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+str(mapping[rname][qname][k][1][1])+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][3])+"\n")	
						elif mapping[rname][qname][i][0][0]< mapping[rname][qname][k][0][0]<=mapping[rname][qname][i][0][1] and mapping[rname][qname][i][0][1]<mapping[rname][qname][k][0][1] and mapping[rname][qname][i][1][1]<mapping[rname][qname][k][1][0]:#ith record rstart > i+1th record rstart
							#rname/pointer_Start/pointer_end/qname/IES_start/IESend
							MDS_IES.write(rname+"\t"+qname+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+str(mapping[rname][qname][i][1][1])+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][i][3])+"\n")		
						else: #rname/mds1_start/mds1_end/mds2_start/mds2_end/qname/MICreadmds1_start.../strand
							if mapping[rname][qname][i][0][0]<=mapping[rname][qname][k][0][0]:
								if not mapping[rname][qname][k][0][1]<=mapping[rname][qname][i][0][1]:
									scrambling.write(rname+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+qname+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][1][1])+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][k][1][1])+"\tplus\n")
							else:
								if not mapping[rname][qname][i][0][1]<=mapping[rname][qname][k][0][1]:
									scrambling.write(rname+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+qname+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][k][1][1])+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][1][1])+"\tplus\n")

					elif mapping[rname][qname][k][3]==mapping[rname][qname][i][3] and mapping[rname][qname][k][3]==0: # same minus strand
						if mapping[rname][qname][k][0][0]< mapping[rname][qname][i][0][0]<=mapping[rname][qname][k][0][1] and mapping[rname][qname][k][0][1]<mapping[rname][qname][i][0][1] and mapping[rname][qname][k][1][0]>mapping[rname][qname][i][1][1]:#ith record rstart > i+1th record rstart
							#rname/pointer_Start/pointer_end/qname/IES_start/IESend
							MDS_IES.write(rname+"\t"+qname+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+str(mapping[rname][qname][i][1][1])+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][i][3])+"\n")	
						elif mapping[rname][qname][i][0][0]< mapping[rname][qname][k][0][0]<=mapping[rname][qname][i][0][1] and mapping[rname][qname][i][0][1]<mapping[rname][qname][k][0][1] and mapping[rname][qname][i][1][0]>mapping[rname][qname][k][1][1]:#ith record rstart > i+1th record rstart
							#rname/pointer_Start/pointer_end/qname/IES_start/IESend
							MDS_IES.write(rname+"\t"+qname+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+str(mapping[rname][qname][k][1][1])+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][3])+"\n")
						else: #rname/mds1_start/mds1_end/mds2_start/mds2_end/qname/MICreadmds1_start.../strand
							if mapping[rname][qname][i][0][0]<=mapping[rname][qname][k][0][0]:
								if not mapping[rname][qname][k][0][1]<=mapping[rname][qname][i][0][1]:
									scrambling.write(rname+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+qname+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][1][1])+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][k][1][1])+"\tminus\n")
							else:
								if not mapping[rname][qname][i][0][1]<=mapping[rname][qname][k][0][1]:
									scrambling.write(rname+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+qname+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][k][1][1])+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][1][1])+"\tminus\n")
					else:
						#rname/mds1_start/mds1_end/mds2_start/mds2_end/qname/MICreadmds1_start.../strand
						if mapping[rname][qname][i][0][0]<=mapping[rname][qname][k][0][0]:
							scrambling.write(rname+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+qname+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][1][1])+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][k][1][1])+"\t"+str(mapping[rname][qname][i][3])+"\t"+str(mapping[rname][qname][k][3])+"\n")
						else:
							scrambling.write(rname+"\t"+str(mapping[rname][qname][k][0][0])+"\t"+str(mapping[rname][qname][k][0][1])+"\t"+str(mapping[rname][qname][i][0][0])+"\t"+str(mapping[rname][qname][i][0][1])+"\t"+qname+"\t"+str(mapping[rname][qname][k][1][0])+"\t"+str(mapping[rname][qname][k][1][1])+"\t"+str(mapping[rname][qname][i][1][0])+"\t"+str(mapping[rname][qname][i][1][1])+"\t"+str(mapping[rname][qname][k][3])+"\t"+str(mapping[rname][qname][i][3])+"\n")

					k=k+1
				i=i+1


	

MDS_IES.close()
scrambling.close()
telocontig.close()
MDS_3_blocks.close()





