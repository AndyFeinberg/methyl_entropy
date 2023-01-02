#!/usr/bin/env python3
import sys
import os
import argparse
import pathlib
import subprocess
import shlex
import xml.etree.ElementTree as ET
def main():
	#Make argparse to get input flags and inputs
	parser = argparse.ArgumentParser(description='This program get a list of files and put them into Bismark pipelines')
	#parser.add_argument("-genome",dest = "gen", help="reference genome selection",default="hg19",required=True)
	parser.add_argument("-fastq_path",dest = "fad", help="path to fastq files",required=True)
	parser.add_argument("-ref_fastq",dest = "fq", help="fastq files",required=True)
	parser.add_argument("-project_dir",dest = "pd", help="project directory",required=True)
	parser.add_argument("-folder",dest = "fd", help="select if make folders",default="T",required=False)
	parser.add_argument("-aligner",dest = "align", help="choose aligner",default="Arioc",required=False)
	parser.add_argument("-seqtype",dest = "type", help="choose sequencing type",default="WGBS",required=False)
	parser.add_argument("-ref_gene",dest = "ref",help="reference genome pathway",required=True)
	parser.add_argument("-qname",dest="qn",help="QNAME section for Arioc",required=False,default="(*:*:*:*:*:*) *")
	#parser.add_argument("-reads",dest="read",help="reads format, single or paired",default="single",required=False)
	args = parser.parse_args()
	#input samples
	samples = args.fq.split(",")
	fadp = spa(args.fad)
	for lf in samples: #loop for each sample

		#reads=lf.split("_")[-1]
		ref_s =lf.split("_")[0]
		#Initialization
		ali_ids= ""
		#Read file in the sample folder
		rf = "ls %s%s*fastq.gz" %(fadp,spa(lf))
		gzff = read_file(rf)
		#check 2 file situation
		gzff_sp = [i.rstrip(".fastq.gz").split("/")[-1].split("_")[0] for i in gzff]
		gzff_uq=list(set(gzff_sp))
		#print(gzff_uq)
		#Sample names
		sample_name=lf.rstrip("/")

		#Directories
		pdp = spa(args.pd)
		sample_dir= pdp+"%s/" %sample_name
		ali_dir = sample_dir+"Align_output/"
		ali_temp_dir=sample_dir+"Align_temp/"
		trim_dir =sample_dir+"Trim_galore_output/"
		output_dir= sample_dir + "bam_output/"
		sub_dir = sample_dir+"sub/"
		#Make directory for each sample
		if args.fd =="T":
			cf(ali_dir)
			cf(ali_temp_dir)
			cf(trim_dir)
			cf(output_dir)
			cf(sub_dir)

		#Loop over each file#
		for fqgz_r in gzff_uq:
			#reconstruct whole fine name use subprocess
			rff = "ls %s%s%s*fastq.gz" %(fadp,spa(lf),fqgz_r)
			fqgz = read_file(rff)
			if len(fqgz)==1:
				reads="single"
				fqgz=fqgz[0]
			elif len(fqgz)==2:
				reads="paired"
				fqgz=fqgz[0]
			else:
				sys.exit("wrong read number")
			#Global filename after Trim galore
			gfn = fqgz.rstrip(".fastq.gz").split("/")[-1]

			#temp folder for alignment
			temp_fd = ali_temp_dir+gfn+"/"

			if args.fd =="T":
				cf(temp_fd)
			#Start write command#
			trim = sub_dir+gfn+"_trim"+".sh"
			#Trim & Move#
			tm=open(trim,"w")
			tm.write(heads(gfn,sample_dir,"10:00:00","24","parallel,shared,skylake","trim"))
			tm.write("%s" %tg(fqgz,trim_dir,reads))
			tm.close()
			tmid = submitp1(trim)
			#Alignmentn
			ali_sub=sub_dir+gfn+"_align"+".sh"
			ali_s=open(ali_sub,"w")


			if args.align == "Bismark":
				ali_s.write(heads(gfn,sample_dir,"26:00:00","24","shared,parallel","BSM_align"))
				#Chagne rgene
				rgene = args.ref.rstrip('/')+'/'
			elif args.align =="Arioc":
				ali_s.write(heads(gfn,sample_dir,"12:00:00","24","gpuk80,gpup100","Arioc_align"))
				rgene =args.ref.rstrip('/')+"/"+ref_s+"_enc"
				ali_s.write("#SBATCH --gres=gpu:4\n")
			if args.align == "Bismark":
				ali_s.write("%s &&" %bsm(ali_dir,gfn,temp_fd,rgene,trim_dir,reads))
			elif args.align =="Arioc":
			##Arioc#
				ali_s.write("%s &&" %aroc(ali_dir,gfn,temp_fd,trim_dir,rgene,reads,args.type,args.qn))
			elif args.align =="RNA":
				ali_s.write(heads(gfn,sample_dir,"6:00:00","24","shared,parallel","RNA_align"))
				ali_s.write("%s &&"%RNA(ali_dir,gfn,temp_fd,trim_dir,rgene,reads))
			#Flagstat# & remove intermidiate sam file
			ali_s.write("%s\n" %(flg(ali_dir,gfn,args.align,args.type)))
			ali_s.close()
			ali_id=submitp2(ali_sub,tmid)
			ali_ids = ali_ids+ali_id+":"

		#Part2 analysis#
		#merge & rmdup
		if reads == "paired":
			rd = "p"
		elif reads=="single":
			rd="s"
		else:
			sys.exit("wrong reads in rmdup")
		part2 = sub_dir+sample_name+"_merge.sh"
		p2 = open(part2,"w")
		if args.type=="WGS":
			par="shared,parallel,lrgmem"
			SNPsplit_out=""
		elif args.type=="WGBS":
			par="lrgmem"
			SNPsplit_out=SNPsplit(reads,ref_s,output_dir,sample_name,args.type,args.ref.rstrip('/'))
		p2.write(heads(sample_name+"_merge",sample_dir,"2-24:00:00","24",par,"merge"))
		p2.write("module load picard\n")
		p2.write("module load gatk/4.0.0\n")

		p2.write("%s " %rmdup(ali_dir,sample_name,output_dir,args.type,rd))
		p2.write("%s \n" %SNPsplit_out)
		p2.close()
		if args.type=="WGBS":
			p3=open(sub_dir+sample_name+"_bsm.sh","w")
			p3.write(heads(sample_name+"_bismark",sample_dir,"1-12:00:00","24","shared,parallel","bsm"))
			p3.write("(bismark_methylation_extractor -%s -o bam_output --gzip --multicore 24 --cytosine_report --no_overlap --bedGraph --genome_folder ~/data/yfang/referenceGenome/hg19 bam_output/%s_cat_dup_rh.bam)&&"%(rd,sample_name))
			p3.write("(rm bam_output/CH*)&&(rm bam_output/CpG*)\n")
			p3.close()
		elif args.type=="WGS":
			#Variant calling
			p3=open(sub_dir+sample_name+"_vcf.sh","w")
			p3.write(heads(sample_name+"_vcf",sample_dir,"24:00:00","24","shared,parallel","bsm"))
			p3.write("(gatk AddOrReplaceReadGroups -I %s%s_cat_dup.bam -O  %s%s_cat_dup_rg.bam -LB lib1 -PU NA -SM H1 -PL illumina)&&" %(output_dir,sample_name,output_dir,sample_name))
			p3.write("(samtools index %s%s_cat_dup_rg.bam') &&" %(output_dir,sample_name))
			p3.write("(gatk HaplotypeCaller -R %s -I %s%s_cat_dup_rg.bam --dbsnp %s/dbsnp.vcf -stand-call-conf 30 -O %s%s.snp.indels.vcf)&&"\
				%(rgene.rstrip("enc")+hg19_Arioc.fa,output_dir,sample_name,args.ref.rstrip('/'),output_dir,sample_name))
			p3.write("(sed 's/MT/M/g' %s%s.snp.index.vcf > %s%s.snp_noheader.vcf) && " %(output_dir,sample_name,output_dir,sample_name))
			p3.write("(vcf_modify2 -vcf %s%s.snp_noheader.vcf -sample %s -output %s) && "%(output_dir,sample_name,sample_name,output_dir))
			#p3.write("(grep -v X %s%s_ASM.vcf >%s%s_ASM_noX.vcf) && (grep -v Y %s%s_ASM_noX.vcf >%s%s_ASM_nosex.vcf) && "%(output_dir,sample_name,output_dir,sample_name,output_dir,sample_name,output_dir,sample_name))
			p3.write("(rm bam_output/CH*)&&(rm bam_output/CpG*)\n")
			p3.close()
		#p2num = submitp2(part2,ali_ids.rstrip(":"))
def RNA(ali,gfn,tmpath,regene,reads):
	if read == "paired":
		R1f = tmpath+gfn+"_val_1.fq.gz"
		R2f = tmpath+gfn.replace("_1","_2")+"_val_2.fq.gz"
		out = "(hisat2 -p 24 --no-softclip --dta -x %s -1 %s -2 %s -S %s/gfn.sam)" %(rgene,R1f,R2f,ali)
	elif read =="single":
		R1f=tmpath+gfn+"_trimmed.fq.gz"
		out ="(hisat2 -p 24 --no-softclip --dta -x %s -U %s -S %s/gfn.sam)" %(rgene,R1f,ali)


	return out

def read_file(cmd):

	check = subprocess.call(cmd,shell=True)
	if check==0:
		gzf = subprocess.check_output(cmd,shell=True)
	else:
		sys.exit("No samples found in fastq directory")
	#Convert files into the list
	out = gzf.decode("utf-8").split("\n")
	out=list(filter(None,out))
	return(out)
def SNPsplit(reads,ref,output_dir,lf,type,ref_dir):
	bamin = "%s%s_cat_dup" %(output_dir,lf)
	out = "(module load gcc/5.5.0) && (module load samtools/1.9)"
	sort_n="(samtools sort -@24 -n -O bam -o %s_n.bam %s.bam) &&" %(bamin,bamin)
	if type=="WGBS":
		if reads=="paired":
			SNPsplit_out="(SNPsplit --paired --no_sort --snp_file %s/out/%s_SNP.txt --bisulfite %s_n.bam)" %(ref_dir,ref,bamin)

		elif reads=="single":
			SNPsplit_out="(SNPsplit --no_sort --snp_file %s/out/%s_SNP.txt --bisulfite %s_n.bam)" %(ref_dir,ref,bamin)

		genome1="&&(samtools sort -@24 -o %s%s_phased.sort.genome1.bam %s_n.genome1.bam) && (samtools index %s%s_phased.sort.genome1.bam)" %(output_dir,lf,bamin,output_dir,lf)
		clean_g1="&& (rm %s_n.genome1.bam)" %bamin
		genome2="&&(samtools sort -@24 -o %s%s_phased.sort.genome2.bam %s_n.genome2.bam) && (samtools index %s%s_phased.sort.genome2.bam)" %(output_dir,lf,bamin,output_dir,lf)
		clean_g2="&& (rm %s_n.genome2.bam)" %bamin
		unassigned="&&(samtools sort -@24 -o %s%s_phased.sort.unassigned.bam %s_n.unassigned.bam) && (samtools index %s%s_phased.sort.unassigned.bam)" %(output_dir,lf,bamin,output_dir,lf)
		clean_ua="&& (rm %s_n.unassigned.bam)" %bamin
		#flagged = "&&(samtools sort -@24 -o %s%s_phased.sort.allele_flagged.bam %s_n.allele_flagged.bam) && (samtools index %s%s_phased.sort.allele_flagged.bam)" %(output_dir,lf,bamin,output_dir,lf)
		clean_fg = "&& (rm %s_n.allele_flagged.bam)" %bamin
		clean_folder="&&(rm -r Align_output/) && (rm Align_temp/*/*sam) && (rm Align_temp/*/*fq) &&"
		SP_out=clean_folder+sort_n+SNPsplit_out+clean_fg+genome1+clean_g1+genome2+clean_g2+unassigned+clean_ua
		#SNPsplit=clean_folder+mbias
	out=out+SP_out
	return(out)
def catmbs(mbias_dir,mbsp):
	out = ("(cat %s*.txt > %s)") %(mbias_dir,mbsp)
	return out
def submitp2(filename,p1ids):
	print("start submit%s" %filename)
	cmd ="sbatch --dependency=afterany:%s %s" %(p1ids,filename)
	status,jobnum=subprocess.getstatusoutput(cmd)
	print("submitting using %s" %cmd)
	if (status ==0):
		print("%s submitted" %filename)
	else:
		sys.exit("Error submitting script")
	return(jobnum.split(" ")[-1])
def submitp1(filename):
	cmd = "sbatch %s" %filename
	status,jobnum = subprocess.getstatusoutput(cmd)
	jobnum = jobnum.split(" ")[-1]
	if (status == 0):
		print("job is %s" % jobnum)
	else:
		sys.exit("Error submitting script")
	return(jobnum)
def spa(spain):
	if spain[-1] == "/":
		out = spain
	else:
		out = spain+"/"
	return(out)

def cf(path): #create folder
	if not os.path.exists(path):
		os.makedirs(path)
	else:
		sys.exit("folder exists")
	return None
def heads(name,sdir,time,ncore,par,jn): #header of the submission script
	out = "#!/bin/sh\n"
	out = out+"#SBATCH --job-name=%s-%s\n" %(jn,name)
	out = out+"#SBATCH --cpus-per-task=%s\n" %ncore
	out = out+"#SBATCH --time=%s\n" %time
	out = out+"#SBATCH -p %s\n" %par
	out = out+"#SBATCH -o %ssub/%s_%s_%%A.out\n" %(sdir,jn,name)
	out = out+"#SBATCH -e %ssub/%s_%s_%%A.err\n" %(sdir,jn,name)
	out = out+"#SBATCH --chdir %s\n"%sdir
	#out=out+ "echo $PWD\n"
	return(out)
def tg(fqgz,tmpath,pair): #Trim Galore
	out = "(module load trim_galore/0.5.0)\n(module load python/3.6-anaconda)\n(module load intel/18.0)\n(module load fastqc)\n"
	#Reconstruct fastq file
	if pair == "paired":
		out = out+"(trim_galore -j 4 --dont_gzip --paired %s %s -o %s)" %(fqgz,fqgz.replace("_1","_2"),tmpath)
		out = out+"&& (fastqc -f fastq -t 24 -o %s %s %s)" %(tmpath,fqgz,fqgz.replace("_1","_2"))
	elif pair =="single":
		out = out+"(trim_galore -j 4 --dont_gzip -o %s %s)" %(tmpath,fqgz)
		
		out=out+"&& (fastqc -f fastq -t 24 -o %s %s)" %(tmpath,fqgz)
	return out
def bsm(ali_dir,gfn,ali_temp_dir,rgene,tmpath,read): #bismark
	if read == "paired":
		R1f = tmpath+gfn+"_val_1.fq.gz"
		R2f = tmpath+gfn.replace("_1","_2")+"_val_2.fq.gz"
		out = "(bismark --bowtie2 --bam --output_dir %s --temp_dir %s %s -1 %s -2 %s)" %(ali_dir,ali_temp_dir,rgene,R1f,R2f)
	elif read =="single":
		R1f=tmpath+gfn+"_trimmed.fq.gz"
		out = "(bismark --bowtie2 --bam --output_dir %s --temp_dir %s %s %s)" %(ali_dir,ali_temp_dir,rgene,R1f)
	return out
def flg(ali_dir,gfn,align,seq): #flagstat
	out = "(module load gcc/5.5.0) && (module load samtools/1.9)&&"
	R1b = ali_dir+gfn+"_%s_%s.bam" %(align,seq)
	out = out+"(samtools flagstat %s)" %R1b
	return(out)
def mbs(mbias_dir,fn):#mbias

	out = "(bismark_methylation_extractor --paired-end --multicore 12 --report --mbias_only --no_overlap --output %s %s)" \
	%(mbias_dir,fn)
	return(out)
def cat(ali_dir,gfn,output_dir):
	out = "(module load gcc/5.5.0) && (module load samtools/1.9)&&"
	mn = "%s%s_cat.bam" %(output_dir,gfn)
	out =out+ "(samtools cat -o %s %s*_val_1_bismark_bt2_pe.bam)" % (mn,ali_dir)
	return(out)
def rmdup(ali_dir,sn,output_dir,align,mode):
	if align == "WGBS":
		mn = "%s*.bam" %(ali_dir)
		out = "(deduplicate_bismark -%s --multiple --bam --output_dir %s %s) && (mv %s*deduplicated.bam %s%s_cat_dup.bam) &&" %(mode,output_dir,mn,output_dir,output_dir,sn)
	elif align =="WGS":
		out = "(module load gcc/5.5.0) && (module load samtools/1.9)&&"
		mg = "%s%s_cat.bam" %(output_dir,sn)#merged
		st="%s%s_cat_st.bam" %(output_dir,sn)
		dup = "%s%s_cat_dup.bam" %(output_dir,sn) #remove duplicated name
		out = out+"(samtools merge -O bam %s %s*_%s.bam) && "%(mg,ali_dir,align)
		out = out + "(samtools sort -@ 22 -o %s %s) && "%(st,mg)
		out=out+"(picard MarkDuplicates I=%s O=%s M=%s%s.txt  REMOVE_DUPLICATES=TRUE TMP_DIR=%s)" %(st,dup,output_dir,sn,output_dir)
		out=out+ "&&(rm %s%s_cat.bam)"%(output_dir,sn)

	return(out)
def ifm(gfn,output_dir,sample_name,sample_dir,p2num,sub_dir,rgene,mpath,mbsp):
	dup = "%s%s_cat_dup.bam" %(output_dir,gfn)
	sort = "%s%s_cat_dup_sort.bam" %(output_dir,gfn)
	out = "(samtools sort -o %s %s) && (samtools index %s) &&" %(sort,dup,sort)
	out = out + "(MARCC_mbias.py -g %s -dup %s -out %s -txt %s -bsmo IFM)\n" %(rgene,dup,mpath,mbsp)
	infom = sub_dir+sample_name+"_IFM.sh"
	infomw = open(infom,"w")
	infomw.write(heads(sample_name+"_IFM",sample_dir,"12:00:00","80G","6","merge_IFM"))
	infomw.write("%s\n" %out)
	infomw.close()
	infomws = submitp2(infom,p2num)
	return(infomws)
def bsmo(gfn,output_dir,rgene,mpath,mbsp,sample_name,sample_dir,p2num,sub_dir): #Bsmooth
	dup = "%s%s_cat_dup.bam" %(output_dir,gfn)
	out = "(MARCC_mbias.py -g %s -dup %s -out %s -txt %s -bsmo bsmooth)" %(rgene,dup,mpath,mbsp)
	bsm = sub_dir+sample_name+"_BSM.sh"
	bsmw = open(bsm,"w")
	bsmw.write(heads(sample_name+"_BSM",sample_dir,"10:00:00","80G","16","merge_BSM"))
	bsmw.write("%s&&(rm methylation_extractor/CH*.txt)&&(rm methylation_extractor/CpG_O*.txt)\n" %out)
	bsmw.close()
	bsmws=submitp2(bsm,p2num)
	return(bsmws)
def aroc(ao_dir,gfn,ao_temp_dir,trim_dir,rgene,read,ali,QN): #Arioc
	#ao_dir: align_output
	#ao_temp_dir=align_temp_dir
	#trim_dir=trim_galore dir
	#ali= alignment type
	out = "(module load gcc/5.5.0)&&(module load python/3.6-anaconda)&&"
	if read == "paired":
		R1f = ao_temp_dir+gfn+"_val_1.fq"
		R1t = trim_dir+gfn+"_val_1.fq"
		R2t = trim_dir+gfn.replace("_1","_2")+"_val_2.fq"
		R2f = ao_temp_dir+gfn.replace("_1","_2")+"_val_2.fq"
		#out = out+"(gunzip %s.gz) &&(gunzip %s.gz)&& " %(R1t,R2t)
		out = out+ "(mv %s %s) && (mv %s %s) && " %(R1t,R1f,R2t,R2f)
		up="P"
	elif read == "single":
		R1f = ao_temp_dir+gfn+"_trimmed.fq"
		R1t = trim_dir+gfn+"_trimmed.fq"
		R2f = "NA"
		#out = out+"(gunzip %s.gz) && (mv %s %s) &&" %(R1t,R1t,R1f)
		#out=""
		up = "U"
	ae=  arocE(gfn,ao_temp_dir,R1f,R2f,read,ali,QN)
	ap = arocP(gfn,ao_temp_dir,rgene,read,ali)

	out = out +"(AriocE %s) && (Arioc%s %s) && (samtools merge -@ 24 -f %s%s_Arioc_uf_%s.bam %s*.sam) &&"\
	%(ae,up,ap,ao_dir,gfn,ali,ao_temp_dir)
	out = out + ("(samtools view -@ 24 -b -F 256 -q 20 -o %s%s_Arioc_%s.bam %s%s_Arioc_uf_%s.bam) && (rm %s%s_Arioc_uf_%s.bam)")\
	%(ao_dir,gfn,ali,ao_dir,gfn,ali,ao_dir,gfn,ali)
	#out=out+"&&(rm -r %s)"%ao_temp_dir #This folder have high disk usage remove after checking
	return out
def arocE(gfn,ao_temp_dir,R1f,R2f,read,ali,QN): #AriocE
	#note for H1 and mouse, used qualityScoreBias=64
	if read == "paired":
		if ali == "WGBS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocE.paired.CT.cfg')
			cfg = ".AE.paired.CT.cfg"
		elif ali == "WGS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocE.paired.cfg')
			cfg=".AE.paired.cfg"
		root=tree.getroot()
		root.find('dataIn').findall('file')[0].text = '%s'%R1f
		root.find('dataIn').findall('file')[1].text = '%s'%R2f
	elif read =="single":
		if ali == "WGBS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocE.unpaired.CT.cfg')
			cfg = ".AE.unpaired.CT.cfg"
		elif ali == "WGS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocE.unpaired.cfg')
			cfg=".AE.unpaired.cfg"
		root=tree.getroot()
		root.find('dataIn').findall('file')[0].text = '%s'%R1f
	root.find('dataIn').set('QNAME',QN)
	root.find('dataOut').find('path').text = '%s'%ao_temp_dir
	out = ao_temp_dir+gfn+cfg
	tree.write(out)
	return out
def arocP(gfn,ao_temp_dir,rgene,read,ali): #AroicP
	# used Vt=G,20,8
	#gfn is the name of outputfile
	#ao_tem_dir is the directory of output file
	#rgene is reference gene name
	if read == "paired":
		if ali == "WGBS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocP.paired.CT.cfg')
			cfg = ".AP.paired.CT.cfg"
		elif ali == "WGS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocP.paired.cfg')
			cfg=".AP.paired.cfg"
		root=tree.getroot()
		root.find('Q').find('paired').findall('file')[0].text=gfn+"_val_1"
		root.find('Q').find('paired').findall('file')[1].text=gfn.replace("_1","_2")+"_val_2"
		root.find('Q').find('paired').set('readID',"*.*")
	elif read =="single":
		if ali == "WGBS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocU.unpaired.CT.cfg')
			cfg = ".AU.unpaired.CT.cfg"
		elif ali == "WGS":
			tree= ET.parse('/work-zfs/afeinbe2/yfang/bin/allele-specific/src/Arioc_temp/AriocU.unpaired.cfg')
			cfg=".AU.unpaired.cfg"
		root=tree.getroot()
		root.find('Q').find('unpaired').findall('file')[0].text=gfn+"_trimmed"
		root.find('Q').find('unpaired').set('readID',"*.*")
	root.find('R').text = '%s'%rgene
	root.find('Q').set("filePath",ao_temp_dir)
	root.find('A').find('sam').text = ao_temp_dir
	out =ao_temp_dir+gfn+cfg
	tree.write(out)
	return out

if __name__=="__main__":main()
