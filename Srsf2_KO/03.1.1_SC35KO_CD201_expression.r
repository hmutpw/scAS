rm(list=ls())
gc()

indir = "../expression_data/download_from_linux/salmon_merge_counts"
#------import gene/isoform information data
gene_infor = read.table("../../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
gene_tpm = read.table(paste(indir,"gene_normalized_TPM.txt",sep="/"),sep="\t",header=T,row.names=1,check.names=F)

dedir = "./diff_res"
test_iso = read.table(paste(dedir,"EC_WT_KO_edgeR_test_res_all.txt",sep="/"),sep="\t")

Procr_test_gene <- test_iso[test_iso$gene_name %in% "Procr", ]
Procr_gene_exp <- as.data.frame(t(gene_tpm["ENSMUSG00000027611",]))
colnames(Procr_gene_exp) <- 'Procr'
Procr_gene_exp$stage <- factor(c("KO","KO","WT","WT","WT"),levels =c("WT","KO"))

plot_frm <- Procr_gene_exp

library(ggplot2)
p = ggplot(plot_frm, aes(x = stage,y = log2(Procr)))+
geom_jitter(shape=16,width=.2)+
ylim(0,8)

p





EC_sig_frm = EC_de_frm[EC_de_frm$PValue<0.05,]
EC_sig_frm_fc = EC_sig_frm[round(abs(EC_sig_frm$logFC),1)>=1,]
EC_sig_frm_fc_high = EC_sig_frm_fc[apply(EC_sig_frm_fc[,3:4],1,function(x){length(x[x>1])>0}),]
write.table(EC_sig_frm_fc_high,paste(outdir,"EC_WT_KO_edgeR_test_sig_res.txt",sep="/"),sep="\t",quote=F)



