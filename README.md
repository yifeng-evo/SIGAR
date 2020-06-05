# SIGAR
SIGAR: **S**plit-read **I**nference of **G**enome **A**rchitecture and **R**earrangements


## Install

Prerequisites of SIGAR:

Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

BWA http://bio-bwa.sourceforge.net

SAMtools http://samtools.sourceforge.net

python https://www.python.org/downloads/

biopython https://biopython.org/wiki/Download

No installation of SIGAR is needed. Just download all the scripts and you are ready to go!
```
git clone https://github.com/yifeng-evo/SIGAR 
```

## Quickstart
```
python SIGAR.py --genome MAC_genome.fasta --reads reads1_u1.fq,reads2_u2.fq  --output_dir SIGAR_output --path_to_SIGAR /path/to/your/sigar/folder/ --contiglist contiglist_MAC_genome.txt --output_name SIGAR -t 12
```

## Manual
1. Check names of paired end reads.

***Be sure to rename paired end reads before running. BWA .sam files automatically remove "/1" "/2" at ends of read names, which will cause a problem for SIGAR to distinguish the two different paired end reads. Quick way to add "_1" or "_2" suffix to read names:***
```
cat reads1.fastq| paste - - - -| sed 's/^\(\S*\)/\1\_1/' | reads1_u1.fastq
cat reads2.fastq| paste - - - -| sed 's/^\(\S*\)/\1\_2/' | reads2_u2.fastq
```

2. Generate a list of interested contigs.

A tab-delimited contig list needs to be provided in the format of "Contig\tcontig_length\n", one contig per line:
e.g.:
```
Contig1	3000
Contig2	5600
...
```
The contig list does not need to include all contigs in the input genome. It should contain at least one line.

3. Run pipeline

Check all parameters by `python SIGAR.py -h`
```
python SIGAR.py --genome MAC_genome.fasta --reads reads1_u1.fq,reads2_u2.fq  --output_dir SIGAR_output --path_to_SIGAR /path/to/your/sigar/folder/ --contiglist contiglist_MAC_genome.txt --output_name SIGAR -t 12
```

4. Filter regions with abnormal coverage (optional, ***highly recommended to remove false discovery***) 

The results of SIGAR on certain genome regions could be noisy due to abnormal high coverage of that region. This could be caused by errors in genome assembly or the high repetitivity of the region.

We are using pileup.sh from BBmap (https://sourceforge.net/projects/bbmap/) to calculate coverage of mapping in BWA files. This script is automatically included when you download SIGAR folder.

```
cd SIGAR_output/bwa
pileup.sh -Xmx6g in=nosecondary_mapq_BWA_MIC_to_MAC.sam out=bwa_cov.txt bincov=bwa_100_cov binsize=100 stdev=t
```
The above command line will output the average converage for every 100bp window in the genome of split reads (specified by `binsize=100`). The input file here nosecondary_mapq_BWA_MIC_to_MAC.sam is the file including all MIC split reads (unmapped reads in bowtie2 end-to-end) mapped to the reference genome. The coverage should follow normal distribution for uniquely-mapped MDS.
This 100bp window can be changed according to the user’s request.

According to the coverage of your genome, you can specify a coverage threshould to filter your results using `coverage_filter.py`.
```
cd SIGAR_output/results
python coverage_filter.py junction_summary_SIGAR_output pointer_summary_SIGAR_output IES_summary_SIGAR_output scrambled_summary_SIGAR_output ../bwa/bwa_100_cov 100 300
```
The first number indicates the binsize you used in pileup.sh. The second number is the threshould you want to use to screen your results. This step will generate a cov_filter file (all results in genome regions below your coverage threshould) and a cov_abnormal file including results that are in noisy regions.

5. Filter results in telomeric regions (optional, ***recommended for low germline DNA coverage dataset*** )

We observed when using whole cell DNA seq or dataset with low germline DNA coverage, a lot of rearrangement junctions/pointers will be inferred in telomeric regions. This could happen for ciliates with flexibility of telomere additions sites. You can filter the results to remove rearrangment features inferred in these regions.


## Interpret Output files

There are 4 tab-delimited output files in the `./SIGAR_output/results/` folder. Their formats are:

```
Junction summary: 1_rname    2_junction_position    3_frequency

Pointer summary: 1_rname    2_pointer_start    3_pointer_end    4_reads# supporting this pointer by 5’ split    5_reads# supporting this pointer by 3’ split    6_pointer_seq    7_comment(IES available or scrambled)

IES summary: 1_rname    2_pointer_start    3_pointer_end    4_pointer_seq    5_IES_seq(without pointer)    6_reads#supporting this IES

Scrambled summary:1_rname    2_MDS1_start    3_MDS1_end    4_MDS2_start    5_MDS2_end    6_query_read_ID    7_query_read_MDS1_start    8_query_read_MDS1_end    9_query_read_MDS2_start    10_query_read_MDS2_end    11_mapping_strand (minus; plus; or 0\t1 means two MDS have different mapping directions)
```

All parameters for the run will be logged in `./SIGAR_output/commands.sh`

The intermediate files are in ./bwa, ./bowtie2 and ./process. The `nosecondary_mapq_BWA_MIC_to_MAC.sam` can be used for reads mapping view using softwares like IGV. 


#08/22/2019:
I haven't tested parse-only function for a while. So I would recommend to start with reads. 
