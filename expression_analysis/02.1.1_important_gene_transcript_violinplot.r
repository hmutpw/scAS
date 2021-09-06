rm(list=ls())
gc()

library(data.table)
geneinfor = read.table("../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T)
row.names(geneinfor) = geneinfor$gene_id
gene_tpm = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",header=T,row.names=1,check.names=F)
#gene_exp_num = apply(gene_tpm,2,function(x){length(x[which(x>1)])})

sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
#stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
stage = c("AEC","HEC","T1_pre_HSC")

samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]
candidates = c("Runx1","Egr1","Spi1","Gata2","Smad4","Gfi1","Nfe2l2","Irf8","Atf2")
candidates = c("Runx1","Egr1","Spi1","Gata2","Smad4","Gfi1","Nfe2l2","Irf8","Atf2")

candidates = c("Sf3b1")

candidate_transcript_infor = transcript_infor[transcript_infor$gene_name %in% candidates,]
#write.table(candidate_transcript_infor,"candidate_transcript_infor.txt",sep="\t",quote=F)
candidate_transcript_tpm = transcript_tpm[intersect(row.names(candidate_transcript_infor), row.names(transcript_tpm)),]
candidate_transcript_tpm_filter = candidate_transcript_tpm[rowMeans(candidate_transcript_tpm)>0,]
candidates_exp_melt = melt(t(candidate_transcript_tpm_filter))
colnames(candidates_exp_melt) = c("sample","transcript_id","TPM")
candidates_exp_melt$stage = sample2stage[as.character(candidates_exp_melt$sample),"stage"]
candidates_exp_melt$transcript_name = candidate_transcript_infor[as.character(candidates_exp_melt$transcript_id),"transcript_name"]


plot_frm = candidates_exp_melt
plot_frm$stage = factor(plot_frm$stage,levels = stage)

line_size = 0.5
text_size = 8
library(ggplot2)
p <- ggplot(plot_frm, mapping = aes(x = stage,y = log2(TPM+1),color = stage))+
geom_violin(position = "dodge", na.rm = TRUE, scale = "width",width = 0.9)+
geom_boxplot(width=0.25,fill="white", size=0.2,outlier.colour=NA, alpha = 0.8)+
geom_jitter(width=0.2,size=0.8,shape=16,color="black", alpha = 0.8)+
labs(x="",y='log2(TPM+1)',title='smart-seq2')+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_wrap(.~transcript_name,nrow=1)
p


candidates = c("Srsf2","Runx1","Gata2","Gfi1")
candidates_geneid = geneinfor[geneinfor$gene_name %in% candidates,]

candidates_exp = gene_tpm[row.names(candidates_geneid),samples]

candidates_exp_melt = melt(t(candidates_exp))
colnames(candidates_exp_melt) = c("sample","gene_id","TPM")

candidates_exp_melt$stage = sample2stage[as.character(candidates_exp_melt$sample),"stage"]
candidates_exp_melt$gene_name = geneinfor[as.character(candidates_exp_melt$gene_id),"gene_name"]
