rm(list=ls())
gc()
library(data.table)
devtools::install_github("ricardo-bion/ggradar",dependencies=TRUE)
library(ggradar)
#####################################################################################
#import events psi data
#####################################################################################
#------psi in a stage less than 5 cells not considered.
event_psi = read.table("../merge_PSI/merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,row.names=1,check.names=F)

#------import events information data
eventInfor = read.table("../merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage[,1]
stage2color = read.table("../../sample_infor/all_color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)
stage = c("AEC","HEC", "T1_pre_HSC","T2_pre_HSC", "E12", "E14", "Adult_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])

#####################################################################################
#plot event type.
#####################################################################################

event_type_stat = t(apply(event_psi,2,function(x){table(eventInfor[names(x[which(!is.na(x))]),"event_type"])}))
event_type_stat_median = apply(event_type_stat,2,function(x){tapply(x,as.character(sample2stage[names(x),"stage"]),median)})
event_type_stat_propor = apply(event_type_stat_median,1,function(x){x/sum(x)})
event_type_melt = melt(t(event_type_stat_propor))
colnames(event_type_melt) = c("stage","event_type","propor")
event_type_melt$event_type = factor(event_type_melt$event_type,levels = c("SE","A3SS","RI","A5SS","MXE"))
event_type_melt$stage = factor(event_type_melt$stage,levels = stage)

plot_frm = event_type_melt

text_size = 8
line_size = 0.5/1.07
library(ggplot2)
p = ggplot(plot_frm, aes(x = "", y = propor, fill = event_type)) + 
geom_bar(stat = "identity", width = 1,color="black") +    
coord_polar(theta = "y") + 
labs(x = "", y = "", title = "")+
#scale_fill_manual(values=light6color)+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major.y = element_line(size=line_size, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_wrap(~stage,ncol=4)

p

event_type_stat_median


