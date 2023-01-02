#!/usr/bin/env python3
import json
import subprocess 
import sys
#Initialize data
print("Starting reformat json")
out_dat = {'STL003':list(),
	'STL011':list(),
	'skin01':list(),
	'skin03':list(),
	'STL001':list(),
	'skin02':list(),
	'H9':list(),
	'HUES64':list(),
	'STL002':list(),
	'HuFGM02':list(),
	'150':list(),
	'149':list(),
	'112':list(),
	'sum':list()}
print("reading in json list")
js_file = open("/work-zfs/afeinbe2/yfang/bin/allele-specific/src/json_file.txt",'r')
file_name=js_file.readlines()
js_file.close()
loc_out = open("SNP_all_xy.txt","w")
hetn=0
for js_name in file_name:
	with open("json/"+js_name.rstrip('\n'), 'r') as js_f:
		js_in = json.load(js_f)
	js_f.close()
	print("Finish reading %s" %js_name)
	
	for js_entry in js_in:
		js_sub = js_entry["AllelicEpigenome-Variant"]['properties']['Subject']['properties']
		if "Summary Metadata" in js_sub.keys():
			ID = js_entry["AllelicEpigenome-Variant"]['value']
			chrom = js_sub["Chromosome"]['value'].lstrip('chr')
			het_CpG = js_sub["Is On Heterogenous CpG"]	
			ref=js_sub["Reference Allele"]['value']
			alt=js_sub["Alternative Allele"]['value']
			pos=js_sub['Position']['value']
			patient_dat = js_sub["Summary Metadata"]["properties"]['Patients With Variant']['items']
			if len(het_CpG)!=0:
				het=het_CpG['value']
				if het =="TRUE":
					hetn=hetn+1
			else:
				het='.'

			loc = '\t'.join([chrom,pos,ID,ref,alt,'100','PASS',"HET="+het,'GT:FI','0/1:1'])
			for pat in patient_dat:
				pat_id=pat['Patient Identifier']['value']
				out_dat[pat_id].append(loc)
			out_dat['sum'].append(loc)
			loc_out.write(loc+"\n")
	
print(hetn)
	

loc_out.close()

print("Finish loading data")

#Write output
temp = open("/work-zfs/afeinbe2/yfang/bin/allele-specific/src/temp_hg19.vcf", "r")
header =  ''.join(temp.readlines())
header2='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
temp.close()
for sample in out_dat:
	file_name=sample+"_xy.vcf"
	print("joining samples")
	
	sample_out='\n'.join(out_dat[sample])
	
	header3=header2+sample+'\n'
	print('writing sample: %s' %sample)
	out_file=open(file_name,"w")
	out_file.write(header+header3+sample_out)
	out_file.close()
	print('Fnish writng sample: %s' %sample)

