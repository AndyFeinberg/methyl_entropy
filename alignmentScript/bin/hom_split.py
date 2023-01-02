#!/usr/bin/env python3
import numpy as np
import argparse
def main():
	parser = argparse.ArgumentParser(description='This program convert all vcf file to the file for ASM pipeline')
	parser.add_argument("-vcf",dest = "vcf", help="vcf file input",required=True)
	parser.add_argument("-sample",dest = "sample", help="sample file in vcf file",required=True)
	parser.add_argument("-output",dest = "opt", help="output directory",required=True)
	args = parser.parse_args()
	sample=args.sample
	vcf_in=open(args.vcf,'r')
	vcf_out=open(args.opt.rstrip('/')+'/'+sample+'_ASM.vcf','w')
	hom_out=open(args.opt.rstrip('/')+'/'+sample+'_hom.vcf','w')
	vcf_in_lines=vcf_in.readlines()
	temp = open("/work-zfs/afeinbe2/yfang/bin/allele-specific/src/temp_hg19.vcf", "r")
	header =  ''.join(temp.readlines())
	header2='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\t'
	temp.close()
	vcf_out.write(header+header2+'\n')
	for lines in vcf_in_lines:
		line_mod=record_modify(lines)
		if line_mod != "NA" | line_mod!="Hom SNP":
			vcf_out.write(line_mod)
		elif line_mod =="Hom SNP":
			hom_out.write(lines)
	vcf_out.close()
	hom_out.close()
def record_modify(record):
	#Make a for loop
	#\t separate: CHROM, POS, ID, REF, ALT, QUAL, FILTER INFO, FORMAT, SAMPLE
	record_sep=record.split('\t')
	#change format section
	fmt=record_sep[8] #format section
	sample=record_sep[9]
	ref=record_sep[3]
	alt=record_sep[4]
	qual=record_sep[5]
	#Check chromosome
	chrom=record_sep[0]
	if chrom.lstrip('chr') != chrom:
		chrom=chrom.lstrip('chr')
	record_sep[0]=chrom
	chr_all=list(range(1,22))
	chr_all=[str(x) for x in chr_all]
	chr_all.extend(['X','Y','MT'])
	if (len(alt) ==1) & (len(ref)==1) & (any([x==chrom for x in chr_all])): # remove all indels and insertions rm all position not in 1:22 chr and XYM
		
		qual=float(qual)
		#Check quality
		if type(qual) is int or type(qual) is float:
			if qual >=10:
				FI="1"
			elif qual<10:
				FI="0"

		GT=sample.split(':')[np.where([x == 'GT' for x in fmt.split(':')])[0][0]]
		
		sample_out=GT+':'+FI
		fmt_out='GT:FI'
		#Check if heterozygous
		if GT[0] != GT[2]:
			#check if phased
			if GT[1] =='|': #phased: find PS section or other in the future?
				PS=sample.split(':')[np.where([x == 'PS' for x in fmt.split(':')])[0][0]]
				sample_out=GT+':'+FI+':'+PS
				fmt_out='GT:FI:PS'
			elif GT[1]=='/': 
				fmt_out='GT:FI'
			else:
				print('Error in process genotype\n')

			#record SNP genotype
			if (GT[0]=="0") & (GT[2] =="1"): #0/1 case or 0|1 case
				geno=ref+'/'+alt
			elif (GT[0] =="1") & (GT[2]=="0"): #1/0 or 1|0
				geno=alt+'/'+ref
			else:
				print('error in record SNP genotype\n')
			geno_out="\t".join([str(i),chrom,record_sep[1],"1",geno])+'\n'
		
			sample_out='0/1'+':'+FI
			fmt_out='GT:FI'
			geno_out='NA'
			record_sep[8]=fmt_out
			record_sep[9]=sample_out
			record_sep[7]='HET=FALSE' #Future function to check if CpG is heterozygous
			record_sep[6]='PASS'
			return("\t".join(record_sep)+"\n")
		elif GT[0] == GT[2]:
			return('Hom SNP')
			
		else:
			print('Error in comparing genotype\n')


	else:
		return('NA')


if __name__=="__main__":main()
