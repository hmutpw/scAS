#!/bin/bash
ref_dir="/mnt/odisk/raw4T/HSC_Autophagy/ref_genome"
#get ERCC transcript fasta.
perl ./03.1.2.get_ERCC_transcript_fasta.pl $ref_dir/ERCC92.gtf $ref_dir/ERCC92.fa $ref_dir/ERCC92.transcript.fa
#merge all wanted fasta files.
cat $ref_dir/gencode.vM22.transcripts.fa $ref_dir/ERCC92.transcript.fa > $ref_dir/gencode.vM22.all.ERCC.transcripts.fa

#cat $ref_dir/gencode.vM22.pc_transcripts.fa $ref_dir/gencode.vM22.lncRNA_transcripts.fa $ref_dir/ERCC92.transcript.fa > $ref_dir/gencode.vM22.pc.lncRNA.ERCC.transcripts.fa
