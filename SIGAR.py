import argparse   ####great! next figure out the 
import sys
parser = argparse.ArgumentParser(description='This program is to infer MDS-MDS junctions and pointers by split-read mapping',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--genome','-g',help='the ciliate MAC genome that you would used for inferring MDS-MDS junctions.',required=True)
parser.add_argument('--reads','-u',help="""input reads, comma-separated file list
paired-end reads are preferred to be merged by programs e.g. bbmap to extend read length
please make sure the name of paired-end reads are NOT SAME, or ending up with /1,/2: BWA sometimes remove /1 and /2 in the output file, which will cause a problem""")
parser.add_argument('--parse_only',action="store_true",help='bowtie2 and bwa results have already been available, you still needs to provide genome and contig list. Also bwa-file and bowtie2-file are required.')
parser.add_argument('--output_dir',default="./",help='the name of output directory')
parser.add_argument('--path_to_sigar',help='the path to this pipeline directory',required=True)
parser.add_argument('--contiglist',help="""the contigs you want to infer junctions/pointers. 
Format: contig name tab contig length per line. e.g.: 
Contig1	10000
Contig2	20000""",required=True)
parser.add_argument('--output_name','-p',default="sigar",help='name of the output files')
parser.add_argument('--bowtie2-mapq',default=3,type=int,help='default: %(default)s,type:%(type)s')
parser.add_argument('--bwa-mapq',default=10,type=int,help='default: %(default)s,type:%(type)s')
parser.add_argument('--bowtie2_file',help='in sam format. It will only be processed when --parse-only is speicified.')
parser.add_argument('--bwa_file',help='in sam format')
parser.add_argument('--clip_length',default=8,type=int,help='the length of hang-over to be considered as a split junction. default: %(default)s,type:%(type)s')
parser.add_argument('--threads','-t',default=1,type=int,help='the number of threads for bowtie2 and bwa running. default: %(default)s,type:%(type)s')
parser.add_argument('--maxfreq',default=100,type=int,help='the pointer file maxmum freq shown. default: %(default)s,type:%(type)s ')
parser.add_argument('--insertion_ths',default=10,type=int,help='The length of (IES+pointer) has to be equal or higher than this threshould. default: %(default)s,type:%(type)s')
parser.add_argument('--mds_align_ths',default=10,type=int,help='The length of mds has to be higher than this number. default: %(default)s,type:%(type)s')
parser.add_argument('--mds_flanking_pointer',default=6,type=int,help='The mds regions have to have least this number of nucleotides aligned other than the pointer seq. default: %(default)s,type:%(type)s')
parser.add_argument('--NM_threshould',default=4,type=int,help='The number of nonmatched nucleotides in bwa alignment. default: %(default)s,type:%(type)s')
parser.add_argument('--mds_rich_ratio',default=0.9,type=float,help='ratio of alignment in a read that a read is considered to be too mds riched and removed from analysis. default: %(default)s,type:%(type)s')
#parser.add_argument('--neglect_teloregion',default=42,help='neglect short ies in these regions')

#deletion_ths=20
#telo_region=36

import os 
import shutil
args=parser.parse_args()
if os.path.isdir(args.output_dir):
	shutil.rmtree(args.output_dir,ignore_errors=True)
os.system("mkdir "+args.output_dir)
f=open(args.output_dir+"/commands.sh","w")
f.write(str(args))
f.close()


os.system("mkdir "+args.output_dir+"/bowtie2/")
os.system("mkdir "+args.output_dir+"/bwa/")
os.system("mkdir "+args.output_dir+"/process/")
os.system("mkdir "+args.output_dir+"/results/")
bwa_mapq=str(args.bwa_mapq)
bowtie2_mapq=str(args.bowtie2_mapq)

if not args.genome:
	print ("MAC genome file is required.")

if not args.parse_only:
	sys.stderr.write("-------------------------bowtie2_database_building------------------------\n")
	os.system("bowtie2-build "+args.genome+" "+args.output_dir+"/bowtie2/genome_bowtie2_index")
	sys.stderr.write("------------------------bowtie2_end_to_end_alignment------------------------\n")
	os.system("bowtie2 --threads "+str(args.threads)+" -x "+args.output_dir+"/bowtie2/genome_bowtie2_index -U "+args.reads+" -S "+args.output_dir+"/bowtie2/bowtie2_end_to_end.sam")
	os.system("samtools view -b -f 4 "+args.output_dir+"/bowtie2/bowtie2_end_to_end.sam | samtools bam2fq - > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq")
	os.system("python "+args.path_to_sigar+"/filter_telo_reads.py "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq "+args.output_dir+"/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq")
	sys.stderr.write("------------------------bwa_database_building------------------------\n")
	os.system("bwa index -p "+args.output_dir+"/bwa/genome_bwa_index "+args.genome)
	sys.stderr.write("------------------------bwa_local_alignment------------------------\n")
	os.system("bwa mem -t "+str(args.threads)+" "+args.output_dir+"/bwa/genome_bwa_index "+args.output_dir+"/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq |samtools view -q "+bwa_mapq+" -h -F 256 - > "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam")
	os.system("samtools view -h -F 4 -q "+bowtie2_mapq+" "+args.output_dir+"/bowtie2/bowtie2_end_to_end.sam > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_mapped_filteredmapq.sam")
else:###what if only bwa/bowtie2 files are provided,especially only bwa, could it still work?
	if args.bowtie2_file:
		os.system("samtools view -b -f 4 "+args.bowtie2_file+" | samtools bam2fq - > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq")
		os.system("samtools view -h -F 4 -q "+bowtie2_mapq+" "+args.bowtie2_file+" > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_mapped_filteredmapq.sam")
		os.system("python "+args.path_to_sigar+"filter_telo_reads.py "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq "+args.output_dir+"/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq")
	else:
		sys.stderr.write("No bowtie2 file provided")
	if args.bwa_file:
		os.system("samtools view -h -q "+bwa_mapq+" -F 256 "+args.bwa_file+" > "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam") 
	else:
		sys.stderr.write("No bwa file provided")
sys.stderr.write("------------------------analysis_starts------------------------\n")
sys.stdout.write("------------------------analysis_starts------------------------\n")
os.system("python "+args.path_to_sigar+"/0_remove_mds_rich_reads.py "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam "+str(args.insertion_ths)+" "+str(args.NM_threshould)+" "+str(args.mds_rich_ratio))
os.system("python "+args.path_to_sigar+"/1_bowtie2_bwa_parse_cigar.py "+args.output_dir+" "+args.contiglist+" "+str(args.clip_length)+" "+str(args.mds_align_ths)+" "+str(args.insertion_ths))
os.system("sort -u -V "+args.output_dir+"/process/junction_output_script1.txt > "+args.output_dir+"/process/sorted_junction_output_script1.txt")
os.system("python "+args.path_to_sigar+"/2_pointer.py "+args.output_dir+" "+str(args.mds_flanking_pointer)+" "+str(args.maxfreq))
os.system("python "+args.path_to_sigar+"/3_chimeric_reads_for_paired_end.py "+args.output_dir+" /bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam_without_mds_rich "+args.contiglist)
os.system("python "+args.path_to_sigar+"/4_chimeric_read_fetch_IES_seq.py "+args.output_dir+" "+args.output_dir+"/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq")
os.system("cat "+args.output_dir+"/process/IES_output_script1.txt "+args.output_dir+"/process/chimeric_ies.txt > "+args.output_dir+"/process/ies_information.txt")
os.system("sort -u -V -t$'\t' "+args.output_dir+"/process/ies_information.txt > "+args.output_dir+"/process/sorted_ies_information.txt")
os.system("python "+args.path_to_sigar+"/5_IES_pipe.py "+args.output_dir+" /results/IES_summary_"+args.output_name+" sorted_ies_information.txt")
os.system("sort -u -V -t$'\t' "+args.output_dir+"/process/IES_output_script1.txt > "+args.output_dir+"/process/sorted_ies_by_insertion_information.txt")
os.system("python "+args.path_to_sigar+"/5_IES_pipe.py "+args.output_dir+" /process/IES_summary_ies_by_insertion sorted_ies_by_insertion_information.txt")
os.system("python "+args.path_to_sigar+"/6_pointer_pipe.py "+args.output_dir+" IES_summary_"+args.output_name+" pointer_summary_"+args.output_name+" "+args.genome)
os.system("python "+args.path_to_sigar+"/7_scrambled.py "+args.output_dir+" "+args.output_name)
os.system("python "+args.path_to_sigar+"/8_pointer_summary_annotate_scramble.py "+args.output_dir+" "+args.output_name)
os.system("python "+args.path_to_sigar+"/9_junction_pipe.py "+args.output_dir+" "+args.output_name)
sys.stderr.write("Thank you for using sigar, check results in the /results/ folder, check bwa and bowtie2 sam files in the /bwa/ and /bowtie2/ folders respectively")
#sys.stdout.write("Thank you for using sigar, check results in the /results/ folder, check bwa and bowtie2 sam files in the /bwa/ and /bowtie2/ folders respectively")
