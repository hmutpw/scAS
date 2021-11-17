#!/bin/bash

for sample in `cat /mnt/odisk/raw4T/HSC_Autophagy/sample_infor/mouse_sample_list.txt`
do
ref_dir="/mnt/odisk/raw4T/HSC_Autophagy/ref_genome"
index_dir="/mnt/odisk/raw4T/HSC_Autophagy/index/STARIndex/mouseGencodeVM22"
fastq_dir="/mnt/hdisk/main4T/tpw_SC_AS/clean_data/$sample"
out_dir="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/alignment_res/$sample"
mkdir -p $out_dir;
echo '[perform read alignment with STAR] '$sample
STAR --runMode alignReads --runThreadN 16 --genomeDir $index_dir --genomeLoad NoSharedMemory --readFilesCommand zcat --sjdbGTFfile $ref_dir/gencode.vM22.primary_assembly.with.ERCC.annotation.gtf --outFileNamePrefix $out_dir/ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --readFilesIn $fastq_dir/$sample".R1.groomed.fq.gz" $fastq_dir/$sample".R2.groomed.fq.gz"
done
