#------import gene/isoform information data
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

stage2color = read.table("../sample_infor/stage2color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)
figure_col = stage2color[stage,"figurecolor"]

gene_level = read.table("gene_level_umap_res.txt",sep="\t",header=T,row.names=1)
isoform_level = read.table("iso_level_umap_res.txt",sep="\t",header=T,row.names=1)

gene_level_tab = data.frame(gene_level,stage = sample2stage[row.names(gene_level),"stage"],type=rep("gene",nrow(gene_level)))
isoform_level_tab = data.frame(isoform_level,stage = sample2stage[row.names(isoform_level),"stage"],type=rep("trancript",nrow(isoform_level)))

plotSNE = rbind(gene_level_tab, isoform_level_tab)
plotSNE$stage = factor(plotSNE$stage, levels = stage)

library(ggplot2)
text_size = 8
line_size = 0.5/1.07

p = ggplot(plotSNE, aes(x = UMAP_1, y = UMAP_2 , color = stage, fill=stage)) + 
geom_point(size = 1,shape=16)+
#xlim(-10,10)+
#ylim(-23,23)+
scale_color_manual(values = figure_col)+
scale_fill_manual(values = figure_col)+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=0.25*line_size, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size,color="black"),
axis.text.x= element_text(size=text_size,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_grid(.~type)

p

