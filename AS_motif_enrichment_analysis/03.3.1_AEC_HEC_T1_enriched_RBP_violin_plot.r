rm(list=ls())
gc()
gene_infor = read.table("../../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
gene_tpm = read.table("../../expression_data/merged_exp/7_stage_gene_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample

stage=c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])


RBP_list = read.table("./rMAPS_enrichment/EHT_enriched_RBPs.txt",sep="\t",header=T)
row.names(RBP_list) = RBP_list$motif_id
RBP_gene_id = unique(as.character(RBP_list$gene_id))

RBP_gene_tpm = gene_tpm[RBP_gene_id,samples]
row.names(RBP_gene_tpm) = gene_infor[row.names(RBP_gene_tpm),"gene_name"]
library(data.table)
RBP_gene_tpm_melt = melt(t(RBP_gene_tpm))
colnames(RBP_gene_tpm_melt) = c("sample","gene_id","TPM")
RBP_gene_tpm_melt$stage = factor(as.character(sample2stage[RBP_gene_tpm_melt$sample,"stage"]),levels = stage)

gene_order = c("Srsf2","Srsf1","Srsf10","Srsf9","Pcbp1","Sfpq","Srsf5","Srsf6","Pcbp2","Srsf3")
RBP_gene_tpm_melt$gene_id = factor(RBP_gene_tpm_melt$gene_id,levels = gene_order)

#------get my theme
text_size = 8
line_size = 0.5/1.07
library(ggplot2)

mytheme = theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA,color=NA),
strip.background = element_rect(fill=NA,color=NA),
strip.text = element_text(size=text_size,color="black"),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
panel.grid.minor = element_line(size=NA, color = NA),
axis.ticks = element_line(size=line_size,color="black"),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')


plot_frm = RBP_gene_tpm_melt

p <- ggplot(plot_frm,aes(x=stage,y=log2(TPM+1), color=stage))+
geom_violin(position = "dodge", na.rm = TRUE, scale = "width",fill=NA,width = 0.9)+
geom_jitter(width=0.25,size= .3,shape=16,color="black")+
scale_y_continuous(limit=c(0,11),breaks=seq(0,10,by=5))+
labs(x="",y='log2(TPM+1)')+
mytheme+
scale_color_manual(values = c("#1e88e5","#c0ca33","#f4511e"))+
facet_wrap(.~gene_id,nrow=2)

p


diff_path = "../../DEAnalysis/diff_res/gene_level/pairwise_test/diff_res"
AEC_HEC_diff_gene = read.table(file.path(diff_path,"AEC_HEC_sig_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_diff_RBP = AEC_HEC_diff_gene[intersect(RBP_gene_id,row.names(AEC_HEC_diff_gene)),]
AEC_HEC_diff_RBP

HEC_T1_diff_gene = read.table(file.path(diff_path,"HEC_T1_pre_HSC_sig_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_diff_RBP = HEC_T1_diff_gene[intersect(RBP_gene_id,row.names(HEC_T1_diff_gene)),]
HEC_T1_diff_RBP

AEC_T1_diff_gene = read.table(file.path(diff_path,"AEC_T1_pre_HSC_sig_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_diff_RBP = AEC_T1_diff_gene[intersect(RBP_gene_id,row.names(AEC_T1_diff_gene)),]
AEC_T1_diff_RBP













