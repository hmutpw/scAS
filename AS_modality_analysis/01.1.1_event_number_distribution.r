rm(list=ls())
gc()
#####################################################################################
#import events psi data
#####################################################################################
#------psi in a stage less than 5 cells not considered.
event_psi = read.table("../merge_PSI/merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,row.names=1,check.names=F)

#------import events information data
eventInfor = read.table("../merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage[,1]
stage2color = read.table("../../sample_infor/stage2color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)
sixStage = c("AEC","HEC", "T1_pre_HSC","T2_pre_HSC", "E12", "E14", "Adult_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% sixStage,])

sixStagecolor = stage2color[sixStage,"figurecolor"]
sixStagefont = stage2color[sixStage,"fontcolor"]
#####################################################################################
#plot all Events PSI distribution during differentiation type1
#####################################################################################
ASENum = apply(event_psi,2,function(x){length(which(!is.na(x)))})

plotFrame = data.frame(sample = names(ASENum),ASENum = ASENum, stage = sample2stage[names(ASENum),'stage'])

sixStagePlotFrame = plotFrame[samples,]
sixStagePlotFrame$stage = factor(sixStagePlotFrame$stage, levels = rev(sixStage))
sortSample = row.names(sixStagePlotFrame[order(sixStagePlotFrame$stage,sixStagePlotFrame$ASENum,decreasing=T),])
sixStagePlotFrame$sample = factor(sixStagePlotFrame$sample,level = sortSample)
sixStagePlotFrame$stage = factor(sixStagePlotFrame$stage, levels = sixStage)
median_num = median(sixStagePlotFrame$ASENum)

#------select wanted stages.
library(ggplot2)
library(gcookbook)
text_size = 8
line_size = 0.5
dodge_width = position_dodge(0.5)
ymin = 0
ymax = 4.1
alpha = 0.25

library(ggplot2)
p <- ggplot(sixStagePlotFrame,aes(x=sample,y=ASENum ,fill=stage))+
geom_bar(stat="identity", color="black", position=dodge_width,size = 0.01)+
geom_hline(yintercept = c(median_num))+
scale_color_manual(values=sixStagefont)+
scale_fill_manual(values=sixStagecolor)+
scale_y_continuous(limits=c(0,9000),breaks=c(0,8000,4000))+
labs(x="",y="number of alternative splicing\nevents in each single cell(X10^3)")+
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
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')


p


















