#!/usr/bin/python
import os, string
from fnmatch import fnmatch
import sys
run_to_check = sys.argv[1]
print "Getting stats for run "+run_to_check
root = '/lustre/scratch115/projects/paxgene/fastQCres_newruns/sorted/'
name_pattern = "fastqc_data.txt"
path_pattern = "*sorted*"
runfiles=[]
filecheck="*"+run_to_check+"*"
for path, subdirs, files in os.walk(root):
    for name in files:
	if fnmatch(path, path_pattern):
        	if fnmatch(name, name_pattern):
			fullpath=os.path.join(path,name)
			if fnmatch(fullpath,filecheck):
				runfiles.append(fullpath)

out_data=[]
for myfile in runfiles:
	with open(myfile,'r') as f:
		for line in f.readlines():
			li=line.lstrip()
			if  li.startswith('#Total Deduplicated'):
				total_dedup_perc= li.split()[3]
				run=myfile.split('/')[7].split('_')[0]
				if run_to_check=="22891" or run_to_check=="21364" or run_to_check=="21121":
					lane=myfile.split('/')[7].split('_')[1].split('.')[0]
					sample=myfile.split('/')[7].split('_')[1].split('.')[1]
				elif run_to_check=="22219" or run_to_check=="22294" or run_to_check=="22226":
                                        lane=myfile.split('/')[7].split('_')[1]
                                        sample=myfile.split('/')[7].split('_')[2].split('.')[0]
				out_data.append ("{}\t{}\t{}\t{}".format(run, lane, sample, total_dedup_perc))			
f.close()
#now save output file
OUTFILE=root+"run" + run_to_check + ".sorted.deduplicated_percentage.txt"
with open(OUTFILE,'w') as outf:
	for item in out_data:
		outf.write("{}\n".format(item,"\n"))
outf.close()
