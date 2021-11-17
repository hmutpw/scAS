#!/bin/bash

for sample in `cat /mnt/odisk/raw4T/HSC_Autophagy/sample_infor/mouse_sample_list.txt`
do
ref_dir="/mnt/odisk/raw4T/HSC_Autophagy/ref_genome"
bam_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/alignment_res/$sample"
out_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/salmon_res/$sample"
mkdir -p $out_dir;
date;
echo "[Convert bam files to sam files]"
samtools view -@ 8 -h $bam_dir/Aligned.toTranscriptome.out.bam > $out_dir/Aligned.toTranscriptome.out.sam
date;
echo "[Filter sam files with transcripts included in fasta files]"
perl 03.2.2.filter_bam_for_salmon.pl $ref_dir/gencode.vM22.all.ERCC.transcripts.fa $out_dir/Aligned.toTranscriptome.out.sam $out_dir/Aligned.toTranscriptome.out.for.salmon.sam
echo '[Perform transcript quantification with Salmon] '$sample
salmon quant -p 8 --libType IU --seqBias --gencode -t $ref_dir/gencode.vM22.all.ERCC.transcripts.fa -a $out_dir/Aligned.toTranscriptome.out.for.salmon.sam -o $out_dir/ 
rm $out_dir/Aligned.toTranscriptome.out.sam
rm $out_dir/Aligned.toTranscriptome.out.for.salmon.sam
done
