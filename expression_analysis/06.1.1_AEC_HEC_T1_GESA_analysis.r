rm(list=ls())
gc()

gene_infor = read.table("../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

gene_tpm_tab = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",header=T,row.names=1,check.names=F)
gene_tpm = gene_tpm_tab[,samples]
gene_tpm_filter = gene_tpm[apply(gene_tpm,1,function(x){length(x[which(x>1)])>=5}),]

gene_tpm_out = data.frame(gene_name = gene_infor[row.names(gene_tpm),"gene_name"],gene_tpm,check.names=F)
gene_tpm_out$gene_name = toupper(gene_tpm_out$gene_name)
write.table(gene_tpm_out,"./GSEA/expression_data/AEC_HEC_T1_tpm.txt",sep="\t",row.names=F,quote=F)

AEC_HEC_T1_stage = t(as.character(sample2stage[samples,"stage"]))
write.table(AEC_HEC_T1_stage,"./GSEA/expression_data/AEC_HEC_T1_tpm.cls",sep="\t",row.names=F,col.names=F,quote=F)

table(AEC_HEC_T1_stage)
length(AEC_HEC_T1_stage)

########################################
#plot GSEA result.
########################################
rm(list=ls())
gc()

splisome_gene = read.table("./GSEA/splisome_important_gene.txt",sep="\t",header=F)
row.names(splisome_gene) = splisome_gene[,1]
colnames(splisome_gene) = c("gene_name","index")

gsea_path = "./GSEA/gsea_out"
AEC_HEC_res = read.table(file.path(gsea_path, "AEC_HEC_KEGG_SPLISOME.txt"),header=T,row.names=1,sep="\t",check.names=F)
HEC_T1_res = read.table(file.path(gsea_path, "HEC_T1_KEGG_SPLISOME.txt"),header=T,row.names=1,sep="\t",check.names=F)
AEC_T1_res = read.table(file.path(gsea_path, "AEC_T1_KEGG_SPLISOME.txt"),header=T,row.names=1,sep="\t",check.names=F)


gsea_res = HEC_T1_res
gsea_tab = gsea_res[,c(1,5,6,7)]
colnames(gsea_tab) = c("gene_name","gene_rank","metric_score","es_score")
gsea_tab$index = as.character(splisome_gene[as.character(gsea_tab$gene_name),"gene_name"])

library(ggplot2)
text_size = 8
line_size = 0.5/1.07
xmax = max(gsea_tab$gene_rank)
ymin = min(gsea_tab$es_score)
ymax = max(gsea_tab$es_score)

p <- ggplot(gsea_tab, mapping = aes(x = gene_rank,y=es_score))+
geom_line(color="green",size=2*line_size)+
scale_x_continuous(limits = c(0,xmax), breaks=seq(0,xmax,by = 10000))+
scale_y_continuous(limits = c(ymin,ymax), breaks=seq(floor(ymin),ceiling(ymax),by = 0.1))+
labs(x="",y='Enrichment score (ES)')+
theme(line = element_line(size=line_size), 
text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p


q <- ggplot(gsea_tab, mapping = aes(x = gene_rank,label=index))+
geom_vline(aes(xintercept=gene_rank,color=index), size = line_size*0.5)+
geom_text(aes(x = gene_rank,y=log10(gene_rank)/10,label=index),na.rm=T)+
scale_x_continuous(limits = c(0,xmax), breaks=seq(0,xmax,by = 10000))+
theme(line = element_line(size=line_size), 
text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
#panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

q

ymin = min(gsea_tab$metric_score)
ymax = max(gsea_tab$metric_score)
r <- ggplot(gsea_tab, mapping = aes(x = gene_rank))+
geom_segment(aes(xend = gene_rank, y = 0,yend = metric_score,color=index))+
geom_text(aes(x = gene_rank,y=log10(gene_rank)/4,label=index),na.rm=T)+
labs(x="",y='Ranked list metric')+
scale_x_continuous(limits = c(0,xmax), breaks=seq(0,xmax,by = 10000))+
scale_y_continuous(limits = c(ymin,ymax), breaks=seq(floor(ymin),ceiling(ymax),by = 0.1))+
theme(line = element_line(size=line_size), 
text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

r

gridExtra::grid.arrange(p,q,r)






