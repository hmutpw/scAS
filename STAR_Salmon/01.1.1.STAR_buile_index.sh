#!/bin/bash
ref_dir="/mnt/odisk/raw4T/HSC_Autophagy/ref_genome"
index_dir="/mnt/odisk/raw4T/HSC_Autophagy/index/STARIndex/mouseGencodeVM22"
STAR --runMode genomeGenerate --runThreadN 16  --sjdbOverhang 100 --genomeDir $index_dir --genomeFastaFiles $ref_dir/GRCm38.primary_assembly.with.ERCC.genome.fa --sjdbGTFfile $ref_dir/gencode.vM22.primary_assembly.with.ERCC.annotation.gtf
