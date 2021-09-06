rm(list=ls())
gc()
#######################################################################
#test for edgeR.
#######################################################################
doedgeR <- function(x, group){
library(edgeR)
x = x
n = ncol(x)/2
group = group
y <- DGEList(counts = x, group = group)
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
if(n == 1){
bcv <- 0.1
ltr <- exactTest(y,dispersion=bcv^2)
}else{
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
ltr <- glmLRT(fit, coef=2)
}
FDR <- p.adjust(ltr$table$PValue, method="BH")
Direction <- decideTestsDGE(ltr,lfc=1)
res<-cbind(ltr$table,FDR,Direction)
res_order = res[order(res[,"FDR"]),]
return(res_order)
}

#######################################################################
#read data from local
#######################################################################
indir = "../expression_data/download_from_linux/salmon_merge_counts"
WT_KO_counts = read.table(paste(indir,"transcript_raw_counts.txt",sep="/"),sep="\t",header=T,row.names=1,check.names=F)
row.names(WT_KO_counts) = substr(row.names(WT_KO_counts),1,18)
#------import gene/isoform information data
transcript_infor = read.table("../../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
transcript_tpm = read.table(paste(indir,"transcript_normalized_TPM.txt",sep="/"),sep="\t",header=T,row.names=1,check.names=F)

#------EC sample
EC_WT_KO_exp = WT_KO_counts[,-3]
EC_condition = factor(c("KO_AGM","KO_AGM","WT_AGM","WT_AGM","WT_AGM"),levels = c("WT_AGM","KO_AGM"))

transcript_median_tpm = t(apply(transcript_tpm[,-3],1,function(x){tapply(x, EC_condition, median)}))
colnames(transcript_median_tpm) = paste(colnames(transcript_median_tpm),"_median_TPM",sep="")

EC_de_res = doedgeR(x = EC_WT_KO_exp, group = EC_condition)
EC_inter_gene = intersect(row.names(EC_de_res),row.names(transcript_infor))
EC_de_frm = data.frame(transcript_infor[EC_inter_gene,],transcript_median_tpm[EC_inter_gene,],EC_de_res[EC_inter_gene,],check.rows=T)
outdir = "./diff_res"
dir.create(outdir, recursive = T)
write.table(EC_de_frm,paste(outdir,"EC_WT_KO_transcript_edgeR_test_res_all.txt",sep="/"),sep="\t",quote=F)
EC_sig_frm = EC_de_frm[EC_de_frm$PValue<0.05,]
EC_sig_frm_fc = EC_sig_frm[round(abs(EC_sig_frm$logFC),1)>=1,]
EC_sig_frm_fc_high = EC_sig_frm_fc[apply(EC_sig_frm_fc[,3:4],1,function(x){length(x[x>1])>0}),]
write.table(EC_sig_frm_fc_high,paste(outdir,"EC_WT_KO_transcript_edgeR_test_sig_res.txt",sep="/"),sep="\t",quote=F)












