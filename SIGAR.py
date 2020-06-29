import argparse   
import sys
parser = argparse.ArgumentParser(description='SIGAR: Split-read Inference of Genome Architecture and Rearrangements',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--genome','-g',help='the ciliate MAC genome (rearranged genome)',required=True)
parser.add_argument('--reads','-u',help="""input reads, comma-separated file list
paired-end reads are recommended to be merged (e.g. by BBmap) to extend read length
please make sure the name of paired-end reads are NOT SAME, or ending up with /1,/2: BWA sometimes removes /1 and /2 in the output file, which will cause a problem""")
parser.add_argument('--parse_only',action="store_true",help='Apply when bowtie2 and bwa results are ready to use, --genome and --contiglist are still required. Use with --bwa_file and --bowtie2_file')
parser.add_argument('--output_dir',default="./",help='the name of output directory')
parser.add_argument('--path_to_sigar',help='the path to this pipeline directory',required=True)
parser.add_argument('--contiglist',help="""A tab-delimited file containing MAC contigs that SIGAR will analyze.
Format: contig name tab contig length. One contig per line. 
e.g.: 
Contig1	2000
Contig2	3000""",required=True)
parser.add_argument('--output_name','-p',default="sigar",help='Suffix of the output files')
parser.add_argument('--bowtie2-mapq',default=3,type=int,help='default: %(default)s,type:%(type)s')
parser.add_argument('--bwa-mapq',default=10,type=int,help='BWA alignments with MAPQ higher than this number are analyzed. default: %(default)s,type:%(type)s')
parser.add_argument('--bowtie2_file',help='in sam format. Use with --parse_only')
parser.add_argument('--bwa_file',help='in sam format. Use with --parse_only')
parser.add_argument('--clip_length',default=8,type=int,help='the minimum length of split read hang-over to be considered. default: %(default)s,type:%(type)s')
parser.add_argument('--threads','-t',default=1,type=int,help='the number of threads for bowtie2 and bwa running. default: %(default)s,type:%(type)s')
parser.add_argument('--maxfreq',default=100,type=int,help='the maxmum number of split reads supporting a pointer on one side that SIGAR will process. It is to accelerate the speed for analyzing MDS_MDS junctions with extremely high coverage. default: %(default)s,type:%(type)s ')
parser.add_argument('--insertion_ths',default=10,type=int,help='The minimum length of (IES+pointer) to be considered. The number should be big enough to distinguish IES from insertions caused by sequencing error or mapping error. default: %(default)s,type:%(type)s')
parser.add_argument('--mds_align_ths',default=10,type=int,help='The minimum length of mds. default: %(default)s,type:%(type)s')
parser.add_argument('--mds_flanking_pointer',default=6,type=int,help='The minimum length of mds regions other than the pointer seq. default: %(default)s,type:%(type)s')
parser.add_argument('--NM_threshold',default=4,type=int,help='The minimum number of nonmatched nucleotides in bwa alignment. default: %(default)s,type:%(type)s')
parser.add_argument('--mds_rich_ratio',default=0.9,type=float,help='A split read labeled with SA tag with MDS ratio higher than this number will be removed from analysis. This is to exclude chimiric reads caused by sequencing error default: %(default)s,type:%(type)s')
parser.add_argument('--BWA_other_args',default="",help='Other BWA parameters separated by space except the default ones (default command: bwa mem -t 1 genome_index reads.fq). e.g. --BWA_other_args="-a". Note BWA is used to align split reads to MAC genome, adding parameters here could change the SIGAR results. ''')
parser.add_argument('--SecondaryFilter',default='yes',help='''BWA alignment files will be filtered by removing secondary alignments, with only primary mapping analyzed.change to 'no' if all alignments need to be analyzed. default: %(default)s''')
parser.add_argument('--Bowtie2_other_args',default="",help='''Other Bowtie2 parameters separated by space except the default ones (--threads,-U,-x,-S) e.g.: --Bowtie2_other_args="-a --fast". Note bowtie2 is only used to remove reads with end-to-end alignment to MAC genome.''')
parser.add_argument('--Telo5',default='C{3,5}A{3,5}C{3,5}A{3,5}',type=str,help='''To remove 5'telomeric MAC reads from bwa input datasets. Similar pattern should be used to specify the telomere string of your MAC genome. The default is to search reads with C(3-5)A(3-5)C(3-5)A(3-5) on 5'end. Use with --TeloRegion. default: %(default)s,type:%(type)s''')
parser.add_argument('--Telo3',default='G{3,5}T{3,5}G{3,5}T{3,5}',type=str,help='''To remove 3'telomeric MAC reads from bwa input datasets. Similar pattern should be used to specify the telomere string of your MAC genome. The default is to search reads with G(3-5)T(3-5)G(3-5)T(3-5) on 3'end. Use with --TeloRegion. default: %(default)s,type:%(type)s''')
parser.add_argument('--TeloRegion',default=25,type=int,help='The region from read ends that SIGAR will search telomere in. Default is searching %(default)sbp sequence from ends,type:%(type)s')
#parser.add_argument('--neglect_teloregion',default=42,help='neglect short ies in these regions')


import os 
import shutil
args=parser.parse_args()
if os.path.isdir(args.output_dir):
	print ("The output directory exists. Will remove the original directory.")
	shutil.rmtree(args.output_dir,ignore_errors=True) # how about false?
os.system("mkdir "+args.output_dir)
f=open(args.output_dir+"/command_line.sh","w")
f.write(str(args))
f.close()


os.system("mkdir "+args.output_dir+"/bowtie2/")
os.system("mkdir "+args.output_dir+"/bwa/")
os.system("mkdir "+args.output_dir+"/process/")
os.system("mkdir "+args.output_dir+"/results/")
bwa_mapq=str(args.bwa_mapq)
bowtie2_mapq=str(args.bowtie2_mapq)


if not args.parse_only:
	sys.stderr.write("-------------------------bowtie2_database_building------------------------\n")
	os.system("bowtie2-build "+args.genome+" "+args.output_dir+"/bowtie2/genome_bowtie2_index")
	sys.stderr.write("------------------------bowtie2_end_to_end_alignment------------------------\n")
	os.system("bowtie2 --threads "+str(args.threads)+" -x "+args.output_dir+"/bowtie2/genome_bowtie2_index -U "+args.reads+" -S "+args.output_dir+"/bowtie2/bowtie2_end_to_end.sam "+args.Bowtie2_other_args)
	os.system("samtools view -b -f 4 "+args.output_dir+"/bowtie2/bowtie2_end_to_end.sam | samtools bam2fq - > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq")
	os.system("python "+args.path_to_sigar+"/filter_telo_reads.py "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq "+args.output_dir+'''/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq "%s" "%s" %d ''' % (args.Telo5,args.Telo3,args.TeloRegion)) # can add info like number of reads
	sys.stderr.write("------------------------bwa_database_building------------------------\n")
	os.system("bwa index -p "+args.output_dir+"/bwa/genome_bwa_index "+args.genome)
	sys.stderr.write("------------------------bwa_local_alignment------------------------\n")
	if args.SecondaryFilter=="yes":
		os.system("bwa mem -t "+str(args.threads)+" "+args.output_dir+"/bwa/genome_bwa_index "+args.output_dir+"/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq "+args.BWA_other_args+" |samtools view -q "+bwa_mapq+" -h -F 256 - > "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam")
	elif args.SecondaryFilter=="no":
		os.system("bwa mem -t "+str(args.threads)+" "+args.output_dir+"/bwa/genome_bwa_index "+args.output_dir+"/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq "+args.BWA_other_args+" |samtools view -q "+bwa_mapq+" -h - > "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam")
	else:
		sys.stderr.write('''SecondaryFilter argument is illegal. Can only be "yes" or "no"''' )
	os.system("samtools view -h -F 4 -q "+bowtie2_mapq+" "+args.output_dir+"/bowtie2/bowtie2_end_to_end.sam > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_mapped_filteredmapq.sam")
else:
	if args.bowtie2_file:
		os.system("samtools view -b -f 4 "+args.bowtie2_file+" | samtools bam2fq - > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq")
		os.system("samtools view -h -F 4 -q "+bowtie2_mapq+" "+args.bowtie2_file+" > "+args.output_dir+"/bowtie2/bowtie2_end_to_end_mapped_filteredmapq.sam")
		os.system("python "+args.path_to_sigar+"filter_telo_reads.py "+args.output_dir+"/bowtie2/bowtie2_end_to_end_unmapped_reads.fq "+args.output_dir+'''/bowtie2/filtered_bowtie2_end_to_end_unmapped_reads.fq "%s" "%s" %d ''' % (args.Telo5,args.Telo3,args.TeloRegion))
	else:
		sys.stderr.write("No bowtie2 file provided")
	if args.bwa_file:
		if args.SecondaryFilter=="yes":
			os.system("samtools view -h -q "+bwa_mapq+" -F 256 "+args.bwa_file+" > "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam") 
		elif args.SecondaryFilter=="no":
			os.system("samtools view -h -q "+bwa_mapq+" "+args.bwa_file+" > "+args.output_dir+"/bwa/nosecondary_mapq_BWA_MIC_to_MAC.sam") 
		else:
			sys.stderr.write('''SecondaryFilter argument is illegal. Can only be "yes" or "no"''' )	
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
sys.stderr.write("Thank you for using SIGAR, check results in the /results/ folder, check bwa and bowtie2 sam files in the /bwa/ and /bowtie2/ folders.")
#sys.stdout.write("Thank you for using SIGAR, check results in the /results/ folder, check bwa and bowtie2 sam files in the /bwa/ and /bowtie2/ folders.")
