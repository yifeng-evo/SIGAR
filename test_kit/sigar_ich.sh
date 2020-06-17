#PBS !/bin/bash
#PBS -l vmem=16G,walltime=24:0:0
#PBS -l nodes=1:ppn=12
module unload python
module load python/2.7.16
module load biopython
module load samtools
module load bowtie2
module load bwa
cd /N/u/yif/Karst/scratch/Ichthyophthirius_weijen/purge
python /N/u/yif/Karst/scratch/SIGAR_github/SIGAR.py --genome uppercase_GCF_000220395.1_JCVI-IMG1-V.1_genomic.fna --reads ./Genomic_Hiseq/cutadapt_G15_1_NoIndex_L001_R1_1.fastq,./Genomic_Hiseq/cutadapt_G15_1_NoIndex_L001_R2_2.fastq  --output_dir /N/u/yif/Karst/scratch/Ichthyophthirius_weijen/sigar --path_to_sigar /N/u/yif/Karst/scratch/SIGAR_github/ --contiglist contigname_ich.txt --output_name ich_sigar -t 12 --Telo5 "C{4,4}A{2,2}C{4,4}A{2,2}C{4,4}A{2,2}" --Telo3 "G{4,4}T{2,2}G{4,4}T{2,2}G{4,4}T{2,2}" --TeloRegion 25 
