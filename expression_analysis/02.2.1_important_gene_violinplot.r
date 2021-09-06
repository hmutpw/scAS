rm(list=ls())
gc()

library(reshape2)
library(data.table)
geneinfor = read.table("../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T)
row.names(geneinfor) = geneinfor$gene_id
gene_tpm = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",header=T,row.names=1,check.names=F)
gene_exp_num = apply(gene_tpm,2,function(x){length(x[which(x>1)])})

sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
stage = c("AEC","HEC","T1_pre_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

stage2color = read.table("../sample_infor/stage2color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)
figure_col = stage2color[stage,"figurecolor"]

##########
#candidates
##########
candidates = c("Pecam1","Cd44","Procr","Runx1","Gfi1","Spn","Ptprc")
candidates = c("Becn1","Sqstm1","Map1lc3a","Map1lc3b","Atg7","Atg12")
ATG_candidates = read.table("D:/HSC_Autophagy/important_genes/ATG_markers.txt",sep="\t",header=T)
row.names(ATG_candidates) = ATG_candidates$gene_id

candidates_id = geneinfor[geneinfor$gene_name %in% candidates,]
candidates_id = row.names(ATG_candidates)


candidate_gene_tpm = gene_tpm[candidates_id,samples]
row.names(candidate_gene_tpm) = geneinfor[row.names(candidate_gene_tpm),"gene_name"]
candidates_exp_melt = melt(t(candidate_gene_tpm))
colnames(candidates_exp_melt) = c("sample","gene_name","TPM")
candidates_exp_melt$stage = sample2stage[as.character(candidates_exp_melt$sample),"stage"]
candidates_exp_melt$gene_name = factor(candidates_exp_melt$gene_name,levels = as.character(ATG_candidates$gene_name))


plot_frm = candidates_exp_melt
plot_frm$stage = factor(plot_frm$stage,levels = stage)

line_size = 0.5
text_size = 8
library(ggplot2)
p <- ggplot(plot_frm, mapping = aes(x = stage,y = log2(TPM+1),color=stage))+
geom_violin(position = "dodge", na.rm = TRUE, scale = "width",width = 0.9,fill=NA)+
#geom_boxplot(width=0.25,fill="white", size=0.2,outlier.colour=NA, alpha = 0.8)+
geom_jitter(width=0.2,size=0.8,shape=16,color="black")+
scale_y_continuous(limits=c(-0.1,10.5),breaks = seq(0,10,by=2.5))+
scale_color_manual(values = figure_col)+
scale_fill_manual(values = figure_col)+
labs(x="",y='log2(TPM+1)',title='important gene expression')+
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
facet_wrap(.~gene_name,nrow=4)
p


candidates = c("Srsf2","Runx1","Gata2","Gfi1")
candidates_geneid = geneinfor[geneinfor$gene_name %in% candidates,]

candidates_exp = gene_tpm[row.names(candidates_geneid),samples]

candidates_exp_melt = melt(t(candidates_exp))
colnames(candidates_exp_melt) = c("sample","gene_id","TPM")

candidates_exp_melt$stage = sample2stage[as.character(candidates_exp_melt$sample),"stage"]
candidates_exp_melt$gene_name = geneinfor[as.character(candidates_exp_melt$gene_id),"gene_name"]
