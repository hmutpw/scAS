
rm(list=ls())
gc()

in_dir = "./AEC_HEC_T1/enrichment_res/AEC_other_HEC_T1_5_stage_all_in_events_enrichment_metascape"
enrichment_res = read.delim(file.path(in_dir,"AEC_other_HEC_T1_5_stage_all_in_events_enrichment_metascape.txt"),sep="\t",header=T,check.names=F)


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
scale_y_continuous(limits=c(0,12),breaks = seq(0,12,by=4))+
mytheme+
coord_flip()
p










