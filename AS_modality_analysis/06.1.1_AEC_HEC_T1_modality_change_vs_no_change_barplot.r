rm(list=ls())
gc()

###########################################################
#get events modality change number.
###########################################################
modality_tab = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1)
stage = c("AEC","HEC","T1_pre_HSC")
stage_mod = modality_tab[,stage]

AEC_HEC_stat = table(stage_mod[,1:2])
AEC_HEC_stat[1,1]=0
HEC_T1_stat = table(stage_mod[,2:3])
HEC_T1_stat[1,1]=0

###########################################################
#get HEC T1 pre_HSC modality change.
###########################################################
modality_tab = AEC_HEC_stat

total_events = sum(modality_tab)
stable_events = sum(modality_tab[row(modality_tab)==col(modality_tab)])
changed_events = total_events-stable_events
plot_frm_1 = data.frame(type=c("change","unchange"),number = c(changed_events, stable_events),stage = "change_unchange")

other_in = sum(modality_tab[,"included"])-modality_tab["included","included"]
other_ex = sum(modality_tab[,"excluded"])-modality_tab["excluded","excluded"]
ex_in_other = sum(modality_tab[c("included","excluded"),])-
sum(modality_tab["included","included"],modality_tab["excluded","excluded"])
other_other = changed_events-c(other_in+other_ex+ex_in_other)

number = c(other_in, other_ex, ex_in_other, other_other)
names(number) = c("other_in", "other_ex","ex_in_other", "other_other")
plot_frm_2 = data.frame(type = names(number), number,stage = rep("EC_T1",4))
plot_frm = rbind(plot_frm_1,plot_frm_2)

text_size = 8
line_size = 0.5
library(ggplot2)
p = ggplot(plot_frm, aes(x = stage, y=number, fill = type))+
geom_bar(stat="identity", width=0.6, colour = "black", size = line_size*0.6)+
scale_y_continuous(limits = c(0,12500),breaks=seq(0,12500,by=2500))+
theme(line = element_line(size=line_size), text = element_text(size=text_size,color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(color = "white",fill=NA),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')
p

plot_frm_rate = tapply(plot_frm$number,plot_frm$stage,function(x){x/sum(x)})
plot_frm_rate

