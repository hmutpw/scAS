#!/bin/bash
for sample in `cat /mnt/odisk/raw4T/HSC_Autophagy/sample_infor/mouse_sample_list.txt`
do
bam_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/alignment_res/$sample"
temp_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/MISO/tempsam/$sample"
mkdir -p $temp_dir
echo "Preparing BAM files Sample: "$sample
samtools view -@ 8 -h $bam_dir/Aligned.sortedByCoord.out.bam > $temp_dir/Aligned.sortedByCoord.out.sam
perl ./01.2.1_samFilter-101bp.pl $temp_dir/Aligned.sortedByCoord.out.sam $temp_dir/Aligned.filter.out.sam
filter_bam_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/MISO/filterbam/$sample"
mkdir -p $filter_bam_dir
samtools view -@ 8 -bS $temp_dir/Aligned.filter.out.sam > $temp_dir/Aligned.filter.out.bam
samtools sort -@ 8 $temp_dir/Aligned.filter.out.bam -o $filter_bam_dir/Aligned.filter.sorted.bam
samtools index -@ 8 $filter_bam_dir/Aligned.filter.sorted.bam
rm -rf $temp_dir;
done
