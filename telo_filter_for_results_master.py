import sys
contig_list=open(sys.argv[1],"r")
pointer=open(sys.argv[2],"r")
IES=open(sys.argv[3],"r")
scrambled=open(sys.argv[4],"r")

pointer_no_telo=open("telo_filter_"+sys.argv[2],"w")
IES_no_telo=open("telo_filter_"+sys.argv[3],"w")
scrambled_no_telo=open("telo_filter_"+sys.argv[4],"w")


telomere_region=50
contig_dict={}
for line in contig_list:
	c_id=line.split("\t")[0]
	c_len=int(line.split("\t")[1])
	contig_dict[c_id]=c_len


def telo_filter(left,right,length,):
	if left>telomere_region and right < length-telomere_region:
		return 1
	else:
		return 0


for line in pointer:
	rname=line.split("\t")[0]
	pointer_start=int(line.split("\t")[1])
	pointer_end=int(line.split("\t")[2])
	pointer_len=pointer_end-pointer_start
	if rname in contig_dict and pointer_len<=20:
		if telo_filter(pointer_start,pointer_end,contig_dict[rname])==1:
			pointer_no_telo.write(line)

for line in IES:
	rname=line.split("\t")[0]
	pointer_start=int(line.split("\t")[1])
	pointer_end=int(line.split("\t")[2])
	if rname in contig_dict and pointer_len<=20:
		if telo_filter(pointer_start,pointer_end,contig_dict[rname])==1:
			IES_no_telo.write(line)


for line in scrambled:
	rname=line.split("\t")[0]
	pos1=int(line.split("\t")[1])
	pos2=int(line.split("\t")[2])
	pos3=int(line.split("\t")[3])
	pos4=int(line.split("\t")[4])
	pos_left=min(pos1,pos2,pos3,pos4)
	pos_right=max(pos1,pos2,pos3,pos4)
	if rname in contig_dict:
		if telo_filter(pos_left,pos_right,contig_dict[rname])==1:
			scrambled_no_telo.write(line)


contig_list.close()
pointer_no_telo.close()
pointer.close()
IES.close()
IES_no_telo.close()
scrambled_no_telo.close()
scrambled.close()
			



