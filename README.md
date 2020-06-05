# SIGAR
SIGAR: Split-read Inference of Genome Architecture and Rearrangements


# Install
Prerequisites of SIGAR:
Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
BWA http://bio-bwa.sourceforge.net
SAMtools http://samtools.sourceforge.net
python https://www.python.org/downloads/
biopython https://biopython.org/wiki/Download

No installation of SIGAR is needed. Just download all the scripts and you are ready to go!
git clone https://github.com/yifeng-evo/SIGAR 

# Quickstart
python SIGAR.py --genome MAC_genome.fasta --reads reads1.fq,reads2,fq  --output_dir SIGAR_output --path_to_SIGAR /path/to/your/sigar/folder/ --contiglist contiglist_MAC_genome.txt --output_name SIGAR -t 12

# Manual
1. Check names of paired end reads.

Be sure to rename paired end reads before running. BWA samfiles automatically removed "/1" "/2" at ends of read names.
quick way to add "_1" or "_2" suffix to read names:

cat file1.fastq| paste - - - -| sed 's/^\(\S*\)/\1\_1/' | file1_u1.fastq

cat file2.fastq| paste - - - -| sed 's/^\(\S*\)/\1\_2/' | file2_u2.fastq

2. Generate contig list.
A tab-delimited contig list needs to be provided in the format of "Contig\tcontig_length\n", one contig per line:
e.g.:
Contig1 3000
Contig2 5600
...
The contig list does not need to include all contigs in the input genome. It should contain at least one line.

3. Run pipeline

4, filter regions with abnormal coverage (optional) 


Output file format:

Junction summary: 1_rname    2_junction_position    3_frequency

Pointer summary: 1_rname    2_pointer_start    3_pointer_end    4_reads# supporting this pointer by 5’ split    5_reads# supporting this pointer by 3’ split    6_pointer_seq    7_comment(IES available or scrambled)

IES summary: 1_rname    2_pointer_start    3_pointer_end    4_pointer_seq    5_IES_seq(without pointer)    6_reads#supporting this IES

Scrambled summary:1_rname    2_MDS1_start    3_MDS1_end    4_MDS2_start    5_MDS2_end    6_query_read_ID    7_query_read_MDS1_start    8_query_read_MDS1_end    9_query_read_MDS2_start    10_query_read_MDS2_end    11_mapping_strand (minus; plus; or 0\t1 means two MDS have different mapping directions)


python SIGAR.py --help to check all the arguments.


#08/22/2019:
I haven't tested parse-only function for a while. So I would recommend to start with reads. 
