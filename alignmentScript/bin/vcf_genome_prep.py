#!/usr/bin/env python3
import pysam
from pysam import VariantFile
import argparse
import os 
import xml.etree.ElementTree as ET
import subprocess
def main():
	parser = argparse.ArgumentParser(description='This program get a list of files and put them into Bismark pipelines')
	parser.add_argument("-vcf",dest = "vcf", help="vcf file input",required=True)
	parser.add_argument("-ref",dest = "ref", help="Folder to referece genome input",required=True)
	parser.add_argument("-strain",dest = "strain", help="desired strain input",required=True)
	parser.add_argument("-build",dest = "gb", help="genome build name",required=True)
	parser.add_argument("-WGS",dest = "wgs", help="whole genome bam file",required=True)
	args = parser.parse_args()
	cwd = os.getcwd()+"/"
	

	#print(args.strain.split(","))
	sample=args.strain
	#for sample in args.strain.split(","):
	print("making folder %s" %sample)
	#Make output folder 
	fd = sample+"/"
	enc=fd+sample+"_enc/"
	Arioc =fd+"Arioc/"
	sub=fd+"sub/"
	out=fd+"out/"
	tmp=fd+"temp/"
	vcf_phase="out/"+sample+".phased.vcf"
	print(vcf_phase)
	vcf_sort = "temp/"+sample+"_st.vcf"
	ref_vcf = "temp/"+sample+"_N_ref.vcf"
	SNP_text="out/"+sample+"_SNP.txt"
	cf(fd)
	cf(enc)
	cf(Arioc)
	cf(sub)
	cf(out)
	cf(tmp)
	print("writing submission files")
	#format AriocEs
	Ariocs(sample,Arioc+"AriocE_%s.gapped.CT.cfg"%sample,'/work-zfs/afeinbe2/yfang/bin/allele-specific/src/ref_gen_temp/AriocE.gapped.CT.cfg',cwd+fd)
	Ariocs(sample,Arioc+"AriocE_%s.nongapped.CT.cfg"%sample,'/work-zfs/afeinbe2/yfang/bin/allele-specific/src/ref_gen_temp/AriocE.nongapped.CT.cfg',cwd+fd)
	#Write submission scripts
	sub1 = open(sub+sample+"_prep.sh","w")
	sub1.write(heads(sample+"_prep","06:00:00","10","parallel,shared",cwd+fd))
	sub1.write("module load picard\n")
	sub1.write("module load python/3.7\n")
	sub1.write("(picard SortVcf I=%s O=%s)&&" %(args.vcf,vcf_sort))
	sub1.write("(whatshap phase -o %s --tag PS --ignore-read-groups %s %s)&&" %(vcf_phase,vcf_sort,args.wgs))
	sub1.write("(vcf_modify.py -vcf %s -SNP_txt %s -ref_vcf %s -sample %s) &&" %(vcf_phase,SNP_text,ref_vcf,sample))
	sub1.write("(SNPsplit_genome_preparation_json --vcf_file %s --reference_genome %s --genome_build %s --strain %s)&&" %(ref_vcf,args.ref,args.gb,sample))
	sub1.write("(rename .N-masked.fa .fna %s_N-masked/*)&&" %(cwd+fd+sample))
	sub1.write("(rm *.txt.gz)&& (rm -r SNPs_%s)" %sample)
	
	sub1.close()
	sub2 = open(sub+sample+"_Arioc.sh","w")
	sub2.write(heads(sample+"_Arioc","4:00:00","24","gpuk80,gpup100",cwd+fd))
	sub2.write("#SBATCH --gres=gpu:4\n")
	
	sub2.write("AriocE Arioc/AriocE_%s.gapped.CT.cfg\n"%sample)
	sub2.write("AriocE Arioc/AriocE_%s.nongapped.CT.cfg\n"%sample)
	sub2.write("cat %s_N-masked/*.fna > out/%s.masked_hg19.fa\n"%(sample,sample))
	sub2.write("samtools faidx out/%s.masked_hg19.fa\n"%(sample))
	sub2.close()
	id_sub = submitp1(sub+sample+"_prep.sh")
	submitp2(sub+sample+"_Arioc.sh",id_sub)
def Ariocs(sample,out,infile,fd):
	tree= ET.parse(infile)
	root=tree.getroot()
	root.find('dataIn').set("filePath",fd+sample+"_N-masked")
	root.find('dataOut').find('path').text=fd+sample+"_enc"
	tree.write(out)
	return None
def heads(name,time,ncore,par,sdir): #header of the submission script
	out = "#!/bin/sh\n"
	out = out+"#SBATCH --job-name=VCF%s\n" %name
	out = out+"#SBATCH --cpus-per-task=%s\n" %ncore
	out = out+"#SBATCH --time=%s\n" %time
	out = out+"#SBATCH -p %s\n" %par
	out = out+"#SBATCH -o %ssub/%s.out\n" %(sdir,name)
	out = out+"#SBATCH -e %ssub/%s.err\n" %(sdir,name)
	out = out+"#SBATCH --chdir %s\n"%sdir
	#out=out+ "echo $PWD\n"
	return(out)
def cf(path): #create folder
	if not os.path.exists(path):
		os.makedirs(path)
	else:
		print("folder exists, not making folder")
	return None
def submitp2(filename,p1ids):
	print("start submit%s" %filename)
	cmd ="sbatch --dependency=afterany:%s %s" %(p1ids,filename)
	status,jobnum=subprocess.getstatusoutput(cmd)
	print("submitting p2 using %s" %cmd)
	if (status ==0):
		print("%s submitted" %filename)
	else:
		sys.exit("Error submitting Part2")
	return(jobnum.split(" ")[-1])
def submitp1(filename):
	cmd = "sbatch %s" %filename
	status,jobnum = subprocess.getstatusoutput(cmd)
	jobnum = jobnum.split(" ")[-1]
	if (status == 0):
		print("job is %s" % jobnum)
	else:
		sys.exit("Error submitting Part1")
	return(jobnum)
if __name__=="__main__":main()
