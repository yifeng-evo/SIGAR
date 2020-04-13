#after getting IES from chimeric_reads
from Bio import SeqIO
import sys
output_path=sys.argv[1]
f=open(output_path+"/process/mds_ies_bwa_same_contig.txt","r")#"MDS_IES_nosecondary_q10_BWA_MIC_to_MAC_Ewoo_sam_same_contig_revised_format.txt"
reads=open(sys.argv[2],"r")#merged_unmapped_sorted_Ewoo_mic_bowtie2.fastq
chimeric_ies=open(output_path+"/process/chimeric_ies.txt","w")#"MDS_IES_nosecondary_q10_BWA_MIC_to_MAC_Ewoo_sam_chimeric_reads_pointer_ies.txt"
#assembly=open("","r")

reads_for_IES={}
contig_for_pointer={}

def reverse_complimentary(seq):
	new_seq=""
	for i in seq:
		if i =="A":
			new_seq="T"+new_seq
		elif i =="T":
			new_seq="A"+new_seq
		elif i =="G":
			new_seq="C"+new_seq
		else:
			new_seq="G"+new_seq
	return new_seq

for line in f:
	qname=line.split("\t")[1]
#	if qname [-2:]=="_1":
#		qname=qname+"/1"
#	if qname [-2:]=="_2":
#		qname=qname+"/2"
	if qname not in reads_for_IES:
		reads_for_IES[qname]=None
f.close()

for record in SeqIO.parse(reads,"fastq"):
	if record.id in reads_for_IES:
		reads_for_IES[record.id]=record.seq

reads.close()
f=open(output_path+"/process/mds_ies_bwa_same_contig.txt","r")
for line in f:
	qname=line.split("\t")[1]
	rname=line.split("\t")[0]
	pointer_start=int(line.split("\t")[2])
	pointer_end=int(line.split("\t")[3])
	pointer_length=pointer_end-pointer_start
	strand=int(line.split("\t")[6])
	IES_start=int(line.split("\t")[4])
	IES_end=int(line.split("\t")[5])
	if qname in reads_for_IES:
	#if qname [-2:]=="_1":
	#	qname=qname+"/1"
	#if qname [-2:]=="_2":
	#	qname=qname+"/2"
		if strand==1:
			IES=reads_for_IES[qname][IES_start:IES_end]
			pointer=reads_for_IES[qname][IES_end:IES_end+pointer_length]
			strand="+"
		else:
			IES=reverse_complimentary(reads_for_IES[qname][IES_start:IES_end])
			pointer=reverse_complimentary(reads_for_IES[qname][IES_end:IES_end+pointer_length])
			temp=IES_start
			IES_start=len(reads_for_IES[qname])-IES_end
			IES_end=len(reads_for_IES[qname])-temp
			strand="-"

		chimeric_ies.write(str(rname)+"\t"+str(pointer_start)+"\t"+str(pointer_end)+"\t"+str(pointer)+"\t"+str(qname)+"\t"+str(IES_start)+"\t"+str(IES_end)+"\t"+str(IES)+"\t"+str(strand)+"\n")


f.close()

chimeric_ies.close()

