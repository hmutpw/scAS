rm(list=ls())
gc()

enrichment_res = read.delim("./AEC_HEC_T1_multi_iso/enrichment_res/775_HEC_multi_isoform_gene_enrichment_filter.txt",sep="\t",header=T,check.names=F)
enrichment_res = read.delim("./AEC_HEC_T1_multi_iso/enrichment_res/1651_HEC_T1_multi_isoform_gene_enrichment_res_filter.txt",sep="\t",header=T,check.names=F)
enrichment_res = read.delim("./AEC_HEC_T1_multi_iso/enrichment_res/1368_T1_only_multi_isoform_gene_enrichment_res_filter.txt",sep="\t",header=T,check.names=F)


enrichment_tab = enrichment_res[,c("Description","LogP")]
enrichment_tab = enrichment_tab[order(enrichment_tab$LogP),]
enrichment_tab$Description = factor(enrichment_tab$Description,levels = rev(as.character(enrichment_tab$Description)))

plot_frm = enrichment_tab
text_size = 8
line_size = 0.5/1.07
library(ggplot2)
#------set plot theme
mytheme = theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA,color=NA),
strip.background = element_rect(fill=NA,color=NA),
strip.text = element_text(size=text_size,color="black"),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
panel.grid.minor = element_line(size=NA, color = NA),
axis.ticks = element_line(size=line_size,color="black"),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, hjust = 0.5,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p <- ggplot(plot_frm,mapping = aes(x=Description,y=-LogP))+
geom_bar(stat="identity",width=0.8,fill="blue")+
mytheme+
coord_flip()
p










