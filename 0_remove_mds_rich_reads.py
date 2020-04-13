import re
import math
import sys

samfile=open(sys.argv[1],"r")
#SAfile=open(sys.argv[1]+"_mds_rich","w")
without_mds_rich=open(sys.argv[1]+"_without_mds_rich","w")
#nm_low=open(sys.argv[1]+"nm_low","w")
insertion_threshold=int(sys.argv[2])
NM_threshold=int(sys.argv[3])
mds_rich_ratio=float(sys.argv[4])
#unmapped_reads=open(sys.argv[2],"r")
#current ratio 0.9
#output=open(sys.argv[3],"w")

def strand(flag):
	plus=1 #plus=1: true; plus=0: reverse. defalut: plus strand
	a=10
	while a >4 and flag >0:
		a=int(math.log(flag,2))
		if a == 4:
			plus=0
		flag=flag-2**a

	return plus

def queryspan(CIGAR, start, strand):
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


	return query_start,query_end,hard_clip_5,hard_clip_3


def NM_parse(NM_section,CIGAR):
	NM=int(NM_section[5:])
	matchM = re.findall(r'(\d+)M', CIGAR)
	intlistM = [int(x) for x in matchM]
	matchI = re.findall(r'(\d+)I', CIGAR)
	intlistI = [int(x) for x in matchI]
	sumI=sum(intlistI)
	if len(intlistI)==len(intlistM)-1 and CIGAR.count("D")==0 and all(k>=insertion_threshold for k in intlistI):
		NM=NM-sumI
	return NM


query_ali={}
query_len={}

for line in samfile:
	if line[0]!="@":
		NM_section=line.split("\t")[11]
		cigar=line.split("\t")[5]
		if NM_parse(NM_section,cigar)<=NM_threshold:
			if "SA:Z:" in line:# to avoid repeats
				qspan=[]
				qname=line.split("\t")[0]
				if qname not in query_ali:
					flag=int(line.split("\t")[1])
					start_pos=int(line.split("\t")[3])
					cigar=line.split("\t")[5]
					seq=line.split("\t")[9]
					qspan=range(queryspan(cigar,start_pos,strand(flag))[0],queryspan(cigar,start_pos,strand(flag))[1])
					SAsection=line.split("\t")[15]#tested with Oxytricha MAC bwa sam file

					for SA in SAsection.split(";"):
						if len(SA)>=3:
							sa_cigar=SA.split(",")[3]
							sa_start=int(SA.split(",")[1])
							if SA.split(",")[2]=="+":
								sa_strand=1
							elif SA.split(",")[2]=="-":
								sa_strand=0
							qspan.extend(range(queryspan(sa_cigar,sa_start,sa_strand)[0],queryspan(sa_cigar,sa_start,sa_strand)[1]))
					
					qspan_len=len(set(qspan))
					query_ali[qname]=qspan_len
					query_length=len(seq)+queryspan(cigar,start_pos,strand(flag))[2]+queryspan(cigar,start_pos,strand(flag))[3]
					query_len[qname]=query_length
				
				ratio=float(query_ali[qname])/float(query_len[qname])
				if ratio<mds_rich_ratio:
					without_mds_rich.write(line)
			else:
				without_mds_rich.write(line)
	else:
		without_mds_rich.write(line)


#SAfile.close()
without_mds_rich.close()
#nm_low.close()

#mds_rich_samfile.close()
#if supplementary region cover the whole read (>90%), throw this read away (notice no matter how low the mapq is )
samfile.close()
#unmapped_reads.close()
