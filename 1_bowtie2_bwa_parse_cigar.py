import re
import math
import sys
output_path=sys.argv[1]
bowtie2=open(output_path+"/bowtie2/bowtie2_end_to_end_mapped_filteredmapq.sam","r")
bwa=open(output_path+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam_without_mds_rich","r") # removed mds-rich reads
telocontig=open(sys.argv[2],"r") #contig list
junction_output=open(output_path+"/process/junction_output_script1.txt","w")
IES_output=open(output_path+"/process/IES_output_script1.txt","w")

clip_length=int(sys.argv[3]) # clip length has to be higher than this threshould
mds_align_ths=int(sys.argv[4])
insertion_ths=int(sys.argv[5])

telolist={}
for line in telocontig:
	telolist[line.split("\t")[0]]=int(line.split("\t")[1])
telocontig.close()

def strand (flag):
	plus="+" #plus=+: true; plus=-: reverse. defalut: plus strand
	a=10
	while a >4 and flag >0:
		a=int(math.log(flag,2))
		if a == 4:
			plus="-"
		flag=flag-2**a

	return plus

def junction_adjust(ref_junction,query_IES_start,query_IES_end,qseq):
	#check if SEQ is reverse complemented.
	pointer_1=0
	pointer_2=0
	if qseq[query_IES_end]==qseq[query_IES_start]: # if pointer on 5' of IES insertion   "TA...."TA
		i=query_IES_end+1
		j=query_IES_start+1
		while j < query_IES_end and i <len(qseq) and qseq[i]==qseq[j]:
			i=i+1
			j=j+1
		pointer_1=qseq[query_IES_end:i]
	elif qseq[query_IES_start-1]==qseq[query_IES_end-1]:
		i=query_IES_end-2
		j=query_IES_start-2
		while i*j>0 and i>=query_IES_start and qseq[i]==qseq[j]:
			i=i-1
			j=j-1
		pointer_2=qseq[i+1:query_IES_end]
	if pointer_1!=0 and pointer_2!=0:
		pointer=pointer_2+pointer_1
		ref_junction=ref_junction-len(pointer_2)
		ref_junction_other_end=ref_junction+len(pointer)
	elif pointer_1==0 and pointer_2!=0:
		pointer=pointer_2
		ref_junction_other_end=ref_junction-len(pointer)
		query_IES_end=query_IES_end-len(pointer) 
	elif pointer_1!=0 and pointer_2==0:
		pointer=pointer_1
		ref_junction_other_end=ref_junction+len(pointer)
		query_IES_start=query_IES_start+len(pointer)
	else:
		pointer=0
		ref_junction_other_end=ref_junction
	IES=qseq[query_IES_start:query_IES_end]

	return min(ref_junction,ref_junction_other_end), max(ref_junction,ref_junction_other_end),pointer, query_IES_start,query_IES_end, IES

def parse_CIGAR(CIGAR,qseq,ref_start):

	matchM = re.findall(r'(\d+)M', CIGAR)
	intlistM = [int(x) for x in matchM]
	sumM=sum(intlistM)
	matchD = re.findall(r'(\d+)D', CIGAR)
	intlistD = [int(x) for x in matchD]
	sumD=sum(intlistD)
	ref_end=sumM+sumD+ref_start
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
	query_end=clip_5+sumM+sumI



	if len(intlistI)==len(intlistM)-1 and CIGAR.count("D")==0 and len(intlistI)>0: # make sure the CIGAR is clean enough to parse
		if all(i>=mds_align_ths for i in intlistM) and all(k>=insertion_ths for k in intlistI):#all MDS and IES have to be higher than the threshold
			ref_junction=[]#junction list
			query_IES_start=[]
			query_IES_end=[]
			ref_junction.insert(0,int(ref_start+intlistM[0]))#initialization
			query_IES_start.insert(0,int(soft_clip_5+intlistM[0]))
			query_IES_end.insert(0,int(query_IES_start[0]+intlistI[0]))
			n=1
			while n<len(intlistI):
				ref_junction.insert(n,int(ref_junction[n-1]+intlistM[n]))
				query_IES_start.insert(n,int(query_IES_end[n-1]+intlistM[n]))
				query_IES_end.insert(n,int(query_IES_start[n]+intlistI[n]))
				n=n+1
			#n=len(intlistI)
			junction_5=[0]*n
			junction_3=[0]*n
			pointer_l=[0]*n
			IES_l=[0]*n
			for z in range(n):
				junction_5[z],junction_3[z],pointer_l[z],query_IES_start[z],query_IES_end[z],IES_l[z]=junction_adjust(ref_junction[z],query_IES_start[z],query_IES_end[z],qseq)
				query_IES_start[z]=query_IES_start[z]+hard_clip_5
				query_IES_end[z]=query_IES_end[z]+hard_clip_5 #query sequence is different from qseq when there's hard clip

			return clip_5,clip_3,junction_5,junction_3,query_IES_start,query_IES_end,pointer_l,IES_l,query_end,ref_end	
		else:
			return clip_5,clip_3,query_end,ref_end
	else:
		return clip_5,clip_3,query_end,ref_end

def determine_L_R_junction(clip_5,clip_3,clip_length,ref_start,ref_end,distance_to_end):
	if clip_3<=clip_length or distance_to_end<=clip_length:
		R_junction="None"
	else:
		R_junction=ref_end
	if clip_5<=clip_length or ref_start<=clip_length:
		L_junction="None"
	else:
		L_junction=ref_start
	return L_junction,R_junction


def parse_samfile(samfile,telolist):
	for line in samfile:
		if line[0]!="@":
			rname=line.split("\t")[2]
			if rname in telolist:
				CIGAR=line.split("\t")[5]
				ref_start=int(line.split("\t")[3])-1 # change to 0_based
				qname=line.split("\t")[0]
				qseq=line.split("\t")[9]
				flag=int(line.split("\t")[1])
				stra=strand(flag)
				if "I" in CIGAR or "S" in CIGAR or "H" in CIGAR:
					parse_result=parse_CIGAR(CIGAR,qseq,ref_start)
					if len(parse_result)==4:# no IES exist:
						clip_5,clip_3,query_end,ref_end=parse_result
						distance_to_end=telolist[rname]-ref_end
						L_junction,R_junction=determine_L_R_junction(clip_5,clip_3,clip_length,ref_start,ref_end,distance_to_end)
						if L_junction!="None" or R_junction!="None":
							junction_output.write(rname+"\t"+qname+"\t"+str(L_junction)+"\t"+str(R_junction)+"\t"+str(ref_start)+"\t"+str(ref_end)+"\t"+str(clip_5)+"\t"+str(query_end)+"\t"+str(stra)+"\n")
					else:
						clip_5,clip_3,junction_5,junction_3,query_IES_start,query_IES_end,pointer_l,IES_l,query_end,ref_end=parse_result
						distance_to_end=telolist[rname]-ref_end
						L_junction,R_junction=determine_L_R_junction(clip_5,clip_3,clip_length,ref_start,ref_end,distance_to_end)
						
						if len(junction_3)==1:
							junction_output.write(rname+"\t"+qname+"\t"+str(L_junction)+"\t"+str(junction_3[0])+"\t"+str(ref_start)+"\t"+str(junction_3[0])+"\t"+str(clip_5)+"\t"+str(query_IES_start[0])+"\t"+str(stra)+"\n")
							#rname/qname/5_JUNCTION/3_JUCNTION/ref_start/ref_end/q_start/q_end/strand
							junction_output.write(rname+"\t"+qname+"\t"+str(junction_5[0])+"\t"+str(R_junction)+"\t"+str(junction_5[0])+"\t"+str(ref_end)+"\t"+str(query_IES_end[0])+"\t"+str(query_end)+"\t"+str(stra)+"\n")
						
						else:
							junction_output.write(rname+"\t"+qname+"\t"+str(L_junction)+"\t"+str(junction_3[0])+"\t"+str(ref_start)+"\t"+str(junction_3[0])+"\t"+str(clip_5)+"\t"+str(query_IES_start[0])+"\t"+str(stra)+"\n")
							i=1
							while i+1<len(junction_5):
								junction_output.write(rname+"\t"+qname+"\t"+str(junction_5[i])+"\t"+str(junction_3[i+1])+"\t"+str(junction_5[i])+"\t"+str(junction_3[i+1])+"\t"+str(query_IES_end[i])+"\t"+str(query_IES_start[i+1])+"\t"+str(stra)+"\n")		
								i=i+1
							junction_output.write(rname+"\t"+qname+"\t"+str(junction_5[i])+"\t"+str(R_junction)+"\t"+str(junction_5[i])+"\t"+str(ref_end)+"\t"+str(query_IES_end[i])+"\t"+str(query_end)+"\t"+str(stra)+"\n")
						
						m=0
						while m<len(junction_5):
							IES_output.write(rname+"\t"+str(junction_5[m])+"\t"+str(junction_3[m])+"\t"+str(pointer_l[m])+"\t"+qname+"\t"+str(query_IES_start[m])+"\t"+str(query_IES_end[m])+"\t"+IES_l[m]+"\t"+str(stra)+"\n")
							m=m+1

parse_samfile(bowtie2,telolist)
parse_samfile(bwa,telolist)
bowtie2.close()
bwa.close()
junction_output.close()
IES_output.close()

