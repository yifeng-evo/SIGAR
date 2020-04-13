import sys

bin_cov=open(sys.argv[5],"r")
junction=open(sys.argv[1],"r")
pointer=open(sys.argv[2],"r")
IES=open(sys.argv[3],"r")
scrambled=open(sys.argv[4],"r")
bin_size=int(sys.argv[6])
max_cov=float(sys.argv[7])
c_junction=open("cov_filter_"+str(max_cov)+sys.argv[1],"w")
c_pointer=open("cov_filter_"+str(max_cov)+sys.argv[2],"w")
c_IES=open("cov_filter_"+str(max_cov)+sys.argv[3],"w")
c_scrambled=open("cov_filter_"+str(max_cov)+sys.argv[4],"w")

junk_junction=open("cov_abn_"+sys.argv[1],"w")
junk_pointer=open("cov_abn_"+sys.argv[2],"w")
junk_IES=open("cov_abn_"+sys.argv[3],"w")
junk_scrambled=open("cov_abn_"+sys.argv[4],"w")


abn_region={}
for line in bin_cov:
	if line[0]!="#":
		coverage=float(line.split("\t")[1])
		if coverage > max_cov:
			rname=line.split("\t")[0]
			position=int(line.split("\t")[2])
			if rname not in abn_region:
				abn_region[rname]=range(position-bin_size,position)
			else:
				for k in range(position-bin_size,position):
					abn_region[rname].append(k)

for line in junction:
	rname=line.split("\t")[0]
	if rname in abn_region:
		jct=int(line.split("\t")[1])
		if jct in abn_region[rname]:
			junk_junction.write(line)
		else:
			c_junction.write(line)
	else:
		c_junction.write(line)

def parse_pointer_file(file,abn_region_list,c_file,junk_file):

	for line in file:
		rname=line.split("\t")[0]

		if rname in abn_region_list:
		#	print abn_region_list[rname]
			pointer_start=int(line.split("\t")[1])
			pointer_end=int(line.split("\t")[2])
		#	print (pointer_end)
		#	print(pointer_start)
			if pointer_start in abn_region_list[rname] or pointer_end in abn_region_list[rname]:
				junk_file.write(line)
			else:
				c_file.write(line)
		else:
			c_file.write(line)

	file.close()
	c_file.close()
	junk_file.close()

parse_pointer_file(pointer,abn_region,c_pointer,junk_pointer)
parse_pointer_file(IES,abn_region,c_IES,junk_IES)


for line in scrambled:
	rname=line.split("\t")[0]
	if rname in abn_region:
		junction_1=int(line.split("\t")[2])
		junction_2=int(line.split("\t")[3])
		junction_3=int(line.split("\t")[4])
		junction_4=int(line.split("\t")[1])
		if junction_1 in abn_region[rname] or junction_2 in abn_region[rname] or junction_3 in abn_region[rname] or junction_4 in abn_region[rname]:
			junk_scrambled.write(line)
		else:
			c_scrambled.write(line)
	else:
		c_scrambled.write(line)

junction.close()
c_junction.close()
junk_junction.close()
scrambled.close()
c_scrambled.close()
junk_scrambled.close()

