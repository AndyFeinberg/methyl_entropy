#!/usr/bin/env python3
import sys
import os
import argparse
import subprocess
import time
def main():
	parser = argparse.ArgumentParser(description='This program get a list of files and download them')
	#parser.add_argument("-genome",dest = "gen", help="reference genome selection",default="hg19",required=True)
	parser.add_argument("-sl",dest = "sample_list", help="path to input lists",required=True)
	parser.add_argument("-d",dest="working_dir",help="path to working directory")
	args = parser.parse_args()
	sample_list=open(args.sample_list,"r")
	file_in=sample_list.read()
	sample_list.close()
	sample_all=file_in.split('\n')
	
	for entry in sample_all:
		print(entry)
		if entry != "":
			element=entry.split('\t')
			
			name=element[0]+'_'+'_'.join(element[1].split(" "))
			reads=element[2]
			SRR=element[3].split(',')
			download_sub(name,reads,SRR,args.working_dir)


def download_sub(name,read,SRR,wd):
	out_dir = name+'_'+read+"/"
	wd=wd.rstrip("/")+"/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	out_dir_temp = out_dir+'out'+'/'
	if not os.path.exists(out_dir_temp):
		os.makedirs(out_dir_temp)
	for fq_list in SRR:
		out_name=out_dir_temp+fq_list+"_dl.sh"
		out1 = open(out_name,"w")
		out1.write(header(out_dir,fq_list,wd))
		sra_temp="/home-4/yfang27@jhu.edu/scratch/yfang/sra_temp/%s.sra*"%fq_list
		out1.write("[ -f %s ] && rm %s\n"%(sra_temp,sra_temp))
		out1.write("prefetch -X 50G -O $PWD %s\n"%fq_list)
		out1.write("fasterq-dump %s.sra -O $PWD -m 80G -t $PWD -e 6 -p -3 --skip-technical -v -P\n"%fq_list)
		out1.write("gzip -1 %s*.fastq\n"%fq_list)
		out1.write("rename %s.sra %s *fastq.gz\n"%(fq_list,fq_list))
		#out1.write("fastq-dump --split-3 --skip-technical -B -W -O $PWD --gzip %s.sra\n" % (fq_list))
		#out1.write("parallel-fastq-dump -t 6 -O $PWD -s %s --tmpdir $PWD --split-3 --skip-technical -B -W --gzip\n" %(fq_list))
		#out1.write("rm %s.sra* \n" % fq_list)
		#out1.write("rm /home-4/yfang27@jhu.edu/scratch/yfang/sra_temp/%s.sra* %s\n"%(fq_list,fq_list))
		out1.close()
		cmd = "sbatch %s" %out_name
		print(cmd)
		time.sleep(1)
		status,jobnum = subprocess.getstatusoutput(cmd)
		if (status == 0):
			print("job is %s" % jobnum)
		else:
			sys.exit("Error submitting Part1")
def header(out_dir,name,wd):
	header=("#!/bin/sh\n")
	header=header+"#SBATCH --job-name=%s_download\n" %name
	header=header+"#SBATCH --cpus-per-task=6\n"
	header=header+"#SBATCH --time=12:00:00\n"
	header=header+"#SBATCH -p shared,parallel,skylake\n"
	header=header+"#SBATCH -o %s%sout/download_%s.out\n" % (wd,out_dir,name)
	header=header+"#SBATCH -e %s%sout/download_%s.err\n"% (wd,out_dir,name)
	header=header+"#SBATCH --chdir %s%s\n"%(wd,out_dir)
	#header=header+"ml python/3.7-anaconda\n"
	return(header)

if __name__=="__main__":main()
