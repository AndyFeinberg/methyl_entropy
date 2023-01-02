
samtools view -s 0.1 -b HUES64_stem_27_undifferentiated_paired_phased.sort.all.bam > HUES64_stem_27_undifferentiated_01_paired_phased.sort.all.bam
samtools view -s 0.3 -b HUES64_stem_27_undifferentiated_paired_phased.sort.all.bam > HUES64_stem_27_undifferentiated_03_paired_phased.sort.all.bam
samtools view -s 0.7 -b HUES64_stem_27_undifferentiated_paired_phased.sort.all.bam > HUES64_stem_27_undifferentiated_07_paired_phased.sort.all.bam
samtools view -s 0.9 -b HUES64_stem_27_undifferentiated_paired_phased.sort.all.bam > HUES64_stem_27_undifferentiated_09_paired_phased.sort.all.bam

samtools view -@20 -s 0.1 -b HUES64_mesoderm_23_paired_phased.sort.all.bam > HUES64_mesoderm_23_01_paired_phased.sort.all.bam
samtools view -@20 -s 0.3 -b HUES64_mesoderm_23_paired_phased.sort.all.bam > HUES64_mesoderm_23_03_paired_phased.sort.all.bam
samtools view -@20 -s 0.5 -b HUES64_mesoderm_23_paired_phased.sort.all.bam > HUES64_mesoderm_23_05_paired_phased.sort.all.bam
samtools view -@20 -s 0.7 -b HUES64_mesoderm_23_paired_phased.sort.all.bam > HUES64_mesoderm_23_07_paired_phased.sort.all.bam
samtools view -@20 -s 0.9 -b HUES64_mesoderm_23_paired_phased.sort.all.bam > HUES64_mesoderm_23_09_paired_phased.sort.all.bam

samtools view -@20 -s 0.1 -b HUES64_endoerm_27_paired_phased.sort.all.bam >HUES64_endoerm_27_01_paired_phased.sort.all.bam
samtools view -@20 -s 0.3 -b HUES64_endoerm_27_paired_phased.sort.all.bam >HUES64_endoerm_27_03_paired_phased.sort.all.bam
samtools view -@20 -s 0.5 -b HUES64_endoerm_27_paired_phased.sort.all.bam >HUES64_endoerm_27_05_paired_phased.sort.all.bam
samtools view -@20 -s 0.7 -b HUES64_endoerm_27_paired_phased.sort.all.bam >HUES64_endoerm_27_07_paired_phased.sort.all.bam
samtools view -@20 -s 0.9 -b HUES64_endoerm_27_paired_phased.sort.all.bam >HUES64_endoerm_27_09_paired_phased.sort.all.bam

samtools merge -@20 -o HUES64_endoerm_01_mesoderm_09_paired_phased.sort.all.bam HUES64_endoerm_27_01_paired_phased.sort.all.bam HUES64_mesoderm_23_09_paired_phased.sort.all.bam
samtools merge -@20 -o HUES64_endoerm_03_mesoderm_07_paired_phased.sort.all.bam HUES64_endoerm_27_03_paired_phased.sort.all.bam HUES64_mesoderm_23_07_paired_phased.sort.all.bam
samtools merge -@20 -o HUES64_endoerm_05_mesoderm_05_paired_phased.sort.all.bam HUES64_endoerm_27_05_paired_phased.sort.all.bam HUES64_mesoderm_23_05_paired_phased.sort.all.bam
samtools merge -@20 -o HUES64_endoerm_07_mesoderm_03_paired_phased.sort.all.bam HUES64_endoerm_27_07_paired_phased.sort.all.bam HUES64_mesoderm_23_03_paired_phased.sort.all.bam
samtools merge -@20 -o HUES64_endoerm_09_mesoderm_01_paired_phased.sort.all.bam HUES64_endoerm_27_09_paired_phased.sort.all.bam HUES64_mesoderm_23_01_paired_phased.sort.all.bam

samtools index -@20  HUES64_endoerm_01_mesoderm_09_paired_phased.sort.all.bam
samtools index -@20  HUES64_endoerm_03_mesoderm_07_paired_phased.sort.all.bam
samtools index -@20  HUES64_endoerm_05_mesoderm_05_paired_phased.sort.all.bam
samtools index -@20  HUES64_endoerm_07_mesoderm_03_paired_phased.sort.all.bam
samtools index -@20  HUES64_endoerm_09_mesoderm_01_paired_phased.sort.all.bam

qsub cpelasm_allele_agnostic.sge HUES64 endoerm_01_mesoderm_09_paired_phased true
qsub cpelasm_allele_agnostic.sge HUES64 endoerm_03_mesoderm_07_paired_phased true
qsub cpelasm_allele_agnostic.sge HUES64 endoerm_05_mesoderm_05_paired_phased true
qsub cpelasm_allele_agnostic.sge HUES64 endoerm_07_mesoderm_03_paired_phased true
qsub cpelasm_allele_agnostic.sge HUES64 endoerm_09_mesoderm_01_paired_phased true


cp ../downstream/output/human_analysis/CPEL_inputs/hg19_all_250bp.gff ~/yfang_dcs04/allele_specific_roadmap_CEPL/work_archive/CpelAsm/data/HUES64/cpelasm/HUES64_allele_agnostic_analysis.gff