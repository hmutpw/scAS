library(data.table)
modality_tab = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1)
modality_stat = apply(modality_tab,2,function(x){table(x[which(x!="")])})
modality = c("included","excluded","middle","bimodal","multimodal")
modality_num = matrix(0,nrow=ncol(modality_tab),ncol = length(modality),
dimnames=list(colnames(modality_tab),modality))
for(i in names(modality_stat)){
for(j in names(modality_stat[[i]])){
modality_num[i,j] = modality_stat[[i]][j]
}
}
sixStages = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
modality_frm = modality_num[sixStages,]
modality_melt_frm = melt(modality_frm)
colnames(modality_melt_frm) = c('stage','modality','number')


plotFrame = modality_melt_frm
plotFrame$stage = factor(plotFrame$stage,levels = sixStages)
plotFrame$modality = factor(plotFrame$modality,levels = modality)
plotFrame$number = as.numeric(plotFrame$number)
colnames(plotFrame) = c('stage','modality','number')

library(plyr)
plot_frm = ddply(plotFrame[,c("stage","modality","number")], .(stage), transform, percent = number/sum(number) )

##########################################################################
#
##########################################################################

text_size = 8
line_size = 0.5
library(ggplot2)
p = ggplot(plot_frm, aes(x = stage, y = percent, fill = modality))+
geom_bar(stat="identity",width=0.7)+
#coord_polar(theta = "y", direction = -1)+
scale_fill_manual(values = c("#c62828","#1565c0","#2e7d32","#6a1b9a","#ef6c00"))+
labs(title = "AS Events Modality Distribution", x = "", y = 'AS Events Number')+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major.y = element_line(size=0.25*line_size, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.text.x= element_text(size=text_size, hjust = 0.5, vjust = 0.5,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')
#+
#facet_wrap(.~stage,nrow=2)

p



plot_frm_in = mean(plot_frm[plot_frm$modality=="included","percent"])
plot_frm_ex = mean(plot_frm[plot_frm$modality=="excluded","percent"])


