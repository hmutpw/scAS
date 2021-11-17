#!/usr/bin/env Rscript
#get salmon result path.
salmon_dir = "/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/salmon_res"
samples = list.files(salmon_dir)
files <- file.path(salmon_dir, samples, "quant.sf")
names(files) <- samples

#get transcript to gene relation.
library(readr)
ref_dir = "/mnt/odisk/raw4T/HSC_Autophagy/ref_genome"
tx2gene <- read_csv(file.path(ref_dir,"gencode.vM22.transcript2gene.with.ERCC.csv"))
head(tx2gene)

#get transcript/gene level counts.
library(tximport)
txi <- tximport(files, type = "salmon", txOut=TRUE)
transcript_counts = txi$counts
transcript_length = txi$length
txi.gene = summarizeToGene(txi, tx2gene = tx2gene)
gene_counts = txi.gene$counts
gene_length = txi.gene$length

#calculate transcript TPM.
enst = grep("^ENS",row.names(transcript_counts),value=T)
spikes = setdiff(row.names(transcript_counts), enst)
transcript_enst_counts = transcript_counts[enst,]
transcript_enst_length = transcript_length[enst,]
transcript_spikes_counts = transcript_counts[spikes,]

library(data.table)
transcript_TPM = data.frame(transcript_id = enst)
row.names(transcript_TPM) = transcript_TPM$transcript_id
for(i in samples){
count_frm = data.frame(counts = transcript_enst_counts[,i],length = transcript_enst_length[,i])
count_norm = apply(count_frm,1,function(x){as.numeric(x[1])/as.numeric(x[2])/1000})
count_tpm = (10^6*count_norm/sum(count_norm))
col_num = ncol(transcript_TPM)
transcript_TPM = data.frame(transcript_TPM, TPM = count_tpm[row.names(transcript_TPM)])
colnames(transcript_TPM)[(col_num+1)] = i
rm(list=c("count_frm","count_norm","count_tpm","col_num"))
}

out_dir = "/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/salmon_merge_counts"
dir.create(out_dir)
write.table(transcript_enst_counts,file.path(out_dir,"transcript_raw_counts.txt"),sep="\t",quote=F)
write.table(transcript_spikes_counts,file.path(out_dir,"transcript_ERCC_counts.txt"),sep="\t",quote=F)
write.table(transcript_TPM,file.path(out_dir,"transcript_normalized_TPM.txt"),sep="\t",quote=F,row.names=F)

#calculate gene TPM.
ensg = grep("^ENS",row.names(gene_counts), value=T)
spikes = setdiff(row.names(gene_counts), ensg)
gene_ensg_counts = gene_counts[ensg,]
gene_ensg_length = gene_length[ensg,]
gene_spikes_counts = gene_counts[spikes,]

gene_TPM = data.frame(gene_id = ensg)
row.names(gene_TPM) = gene_TPM$gene_id
for(i in samples){
count_frm = data.frame(counts = gene_ensg_counts[,i],length = gene_ensg_length[,i])
count_norm = apply(count_frm,1,function(x){as.numeric(x[1])/as.numeric(x[2])/1000})
count_tpm = (10^6*count_norm/sum(count_norm))
col_num = ncol(gene_TPM)
gene_TPM = data.frame(gene_TPM, TPM = count_tpm[row.names(gene_TPM)])
colnames(gene_TPM)[(col_num+1)] = i
rm(list=c("count_frm","count_norm","count_tpm","col_num"))
}

write.table(gene_ensg_counts,file.path(out_dir,"gene_raw_counts.txt"),sep="\t",quote=F)
write.table(gene_TPM,file.path(out_dir,"gene_normalized_TPM.txt"),sep="\t",quote=F,row.names=F)
write.table(gene_spikes_counts,file.path(out_dir,"gene_ERCC_counts.txt"),sep="\t",quote=F)
