#!/bin/bash

for sample in `cat /mnt/odisk/raw4T/HSC_Autophagy/sample_infor/mouse_sample_list.txt`
do
bam_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/MISO/filterbam/$sample"
out_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/MISO/miso_res"
index_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/MISO/preFileForMISO"
readLen=101

echo "Running Sample "$sample
date;
mkdir -p $index_dir/estimateInsertSize/insert-dist/$sample
echo 'estimating insert size and standard deviation for MISO.'
pe_utils --compute-insert-len $bam_dir/Aligned.filter.sorted.bam $index_dir/estimateInsertSize/gencode.vM22.primary_assembly.annotation.min_1000.const_exons.gff --output-dir $index_dir/estimateInsertSize/insert-dist/$sample
mean=`head -1 $index_dir/estimateInsertSize/insert-dist/$sample/Aligned.filter.sorted.bam.insert_len | awk -F '[=|,]' '{print $2}'`
sd=`head -1 $index_dir/estimateInsertSize/insert-dist/$sample/Aligned.filter.sorted.bam.insert_len | awk -F '[=|,]' '{print $4}'`

sample_out=$out_dir/sample_output/$sample
mkdir -p $sample_out
date;
echo "Running MISO with Sample: "$sample
miso --run $index_dir/gffIndex/A3SS $bam_dir/Aligned.filter.sorted.bam --output-dir $sample_out/A3SS/ --read-len $readLen --paired-end $mean $sd --event-type=A3SS
miso --run $index_dir/gffIndex/A5SS $bam_dir/Aligned.filter.sorted.bam --output-dir $sample_out/A5SS/ --read-len $readLen --paired-end $mean $sd --event-type=A5SS
miso --run $index_dir/gffIndex/MXE $bam_dir/Aligned.filter.sorted.bam --output-dir $sample_out/MXE/ --read-len $readLen --paired-end $mean $sd --event-type=MXE
miso --run $index_dir/gffIndex/RI $bam_dir/Aligned.filter.sorted.bam --output-dir $sample_out/RI/ --read-len $readLen --paired-end $mean $sd --event-type=RI
miso --run $index_dir/gffIndex/SE $bam_dir/Aligned.filter.sorted.bam --output-dir $sample_out/SE/ --read-len $readLen --paired-end $mean $sd --event-type=SE

date;
summary_out=$out_dir/summary_output/$sample
mkdir -p $summary_out
echo "Summarizing MISO with Sample: "$sample
summarize_miso --summarize-samples $sample_out/A3SS/ $summary_out/A3SS/
summarize_miso --summarize-samples $sample_out/A5SS/ $summary_out/A5SS/
summarize_miso --summarize-samples $sample_out/MXE/ $summary_out/MXE/
summarize_miso --summarize-samples $sample_out/RI/ $summary_out/RI/
summarize_miso --summarize-samples $sample_out/SE/ $summary_out/SE/
rm -rf $index_dir/estimateInsertSize/insert-dist/$sample
done
