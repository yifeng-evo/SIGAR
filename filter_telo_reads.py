import sys
from Bio import SeqIO
import re

unmapped_fastq=sys.argv[1]
unmapped_fastq_filtered=open(sys.argv[2],"w")
telomere_5=sys.argv[3]
telomere_3=sys.argv[4]
telomere_region=int(sys.argv[5])


CAtelo=r"%s" % (str(telomere_5))
GTtelo=r"%s" % (str(telomere_3))#r"G{3,5}T{3,5}G{3,5}T{3,5}"

k=0
for record in SeqIO.parse(unmapped_fastq,"fastq"):
	if re.search(CAtelo,str(record.seq[:telomere_region])) or re.search(GTtelo,str(record.seq[-telomere_region:])):
		pass
	else:
		k=k+1
		SeqIO.write(record,unmapped_fastq_filtered,"fastq")


print ("Number of reads input for BWA local mapping %d" %(k))
unmapped_fastq_filtered.close()
