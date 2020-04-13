import sys
from Bio import SeqIO
import re

unmapped_fastq=sys.argv[1]
unmapped_fastq_filtered=open(sys.argv[2],"w")

CAtelo=r"C{3,5}A{3,5}C{3,5}A{3,5}"
GTtelo=r"G{3,5}T{3,5}G{3,5}T{3,5}"

for record in SeqIO.parse(unmapped_fastq,"fastq"):
	if re.search(CAtelo,str(record.seq[:25])) or re.search(GTtelo,str(record.seq[-25:])):
		pass
	else:
		SeqIO.write(record,unmapped_fastq_filtered,"fastq")

unmapped_fastq_filtered.close()
