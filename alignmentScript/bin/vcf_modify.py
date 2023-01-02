#!/usr/bin/env python3
import pysam
from pysam import VariantFile
import argparse
import os 
import xml.etree.ElementTree as ET
import subprocess
def main():		
	parser = argparse.ArgumentParser(description='This program get a list of files and put them into vcf file for SNPsplit')
	parser.add_argument("-vcf",dest = "vcf", help="vcf file input",required=True)
	parser.add_argument("-sample",dest = "sample", help="sample name",required=True)
	parser.add_argument("-ref_vcf",dest = "refv", help="reference vcf file",required=True)
	parser.add_argument("-SNP_txt",dest = "snptxt", help="reference vcf file",required=True)
	args = parser.parse_args()
	hd_form=list()
	sample =args.sample
	print("writing vcf and SNP")
	myvcf=VariantFile(args.vcf,"r")
	#Add header
	for formats in myvcf.header.formats:
		hd_form.append(formats)
	#Check if header contain format FI
	if not "FI" in hd_form:
		myvcf.header.formats.add("FI","1","Integer","Whether a sample was a Pass(1) or fail (0) based on FILTER values")
	ref_vcf = args.refv
	SNP_txt=args.snptxt
	vcf_out = VariantFile(ref_vcf,'w', header=myvcf.header)
	SNP_out = open(SNP_txt,"w")
	i=0
	for variant in myvcf:
		#Finding heterozygous SNP
		if variant.samples[sample]['GT'][0] != variant.samples[sample]['GT'][1]:

			geno_comb=list(variant.ref)+list(variant.alts)
			#record real genotype as in the vcf file, 0/1 = ref/alt or 1/0 = alt/ref etc
			geno=geno_comb[variant.samples[sample]['GT'][0]]+"/"+geno_comb[variant.samples[sample]['GT'][1]]
			variant.samples[sample]['GT']=(1,1)
			i=i+1
		else:
			variant.samples[sample]['GT']=(0,1)

		if type(variant.qual) is int or type(variant.qual) is float:
			if variant.qual >=10:
				variant.samples[sample]['FI']=1
			elif variant.qual<10:
				variant.samples[sample]['FI']=0
		#check the format of chromosome
		if variant.chrom.lstrip("chr")==variant.chrom:#without chr
			chrom=variant.chrom
		else:
			chrom=variant.chrom.lstrip("chr")
		if len(geno)==3:		
			SNP_out.write("\t".join([str(i),chrom,str(variant.pos),"1",geno])+"\n")
			vcf_out.write(variant)
	vcf_out.close()
	SNP_out.close()
if __name__=="__main__":main()
