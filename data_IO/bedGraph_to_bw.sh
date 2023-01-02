#!/bin/bash
#$ -l cee,mf=1G,h_vmem=1.5G
#$ -pe local 1
#$ -o "../logfiles/bedGrpah-bw-$JOB_ID.err"
#$ -j y
#$ -m ebas

#$ -N bedGraph-bw
#$ -cwd
# Args
bedGraph_in="$1"
output_dir="$2"
intermediate_dir="$3"
chrom_size_file="$4"
fn=${bedGraph_in/*\//}
#bw_out=${fn/\.bedGraph/.bw}

# Cut bedfile for first row (data)
cut -f1,2,3,4 $bedGraph_in > $intermediate_dir$fn.cut
bedGraphToBigWig $intermediate_dir$fn.cut $chrom_size_file $output_dir$bw_out

#submission
#for fn in ../data_submission/bedGraph_file/compliment_MML_NME_model_mouse/*{mml,nme}*; 
#do echo qsub data_IO/bedGraph_to_bw.sh $fn ${fn//bedGraph/bw} ../data_submission/bedGraph_to_bw_intermediate/ ../data_submission/mm10.chrom.size;
#done