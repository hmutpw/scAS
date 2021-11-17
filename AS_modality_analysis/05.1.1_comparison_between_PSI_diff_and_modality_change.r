rm(list=ls())
gc()
##################################
#------get sample infor.
##################################
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T,check.names=F,stringsAsFactors=F)
row.names(sample2stage) = sample2stage[,1]
psi_path = "../merge_PSI/merged_PSI/"
event_psi = read.table(file.path(psi_path,"all_PSI_Tab.txt"),sep="\t",header=T,row.names=1,check.names=F)

stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])
stage_psi = event_psi[,samples]
stage_psi_cv = t(apply(stage_psi,1,function(x){tapply(x,sample2stage[names(x),"stage"],function(x){sd(x,na.rm=T)/mean(x,na.rm=T)})}))
stage_psi_median = t(apply(stage_psi,1,function(x){tapply(x,sample2stage[names(x),"stage"],function(x){median(x,na.rm=T)})}))

##################################
#------get modality change
##################################
modality_path = "../Anchor/EventModality"
modality_tab = read.table(file.path(modality_path, "7_stage_event_modality.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_mod_tab = modality_tab[,c("AEC","HEC")]
AEC_HEC_mod_filter = AEC_HEC_mod_tab[apply(AEC_HEC_mod_tab,1,function(x){length(x[x==""])<length(x)}),]
AEC_HEC_mod_change = AEC_HEC_mod_filter[apply(AEC_HEC_mod_filter,1,function(x){x[1]!=x[2]}),]
AEC_HEC_mod_change_no_others = AEC_HEC_mod_change[apply(AEC_HEC_mod_change,1,function(x){
length(x[x==""])==1}),]
AEC_HEC_mod_change_others_others = AEC_HEC_mod_change[setdiff(row.names(AEC_HEC_mod_change),
row.names(AEC_HEC_mod_change_no_others)),]
mod_change_type = c(rep("no_others",nrow(AEC_HEC_mod_change_no_others)),
rep("other_other",nrow(AEC_HEC_mod_change_others_others)))
names(mod_change_type) = c(row.names(AEC_HEC_mod_change_no_others),row.names(AEC_HEC_mod_change_others_others))

AEC_HEC_no_mod_change = AEC_HEC_mod_filter[apply(AEC_HEC_mod_filter,1,function(x){x[1]==x[2]}),]
event_change = c(rep("mod_change",nrow(AEC_HEC_mod_change)),rep("no_mod_change",nrow(AEC_HEC_no_mod_change)))
names(event_change) = c(row.names(AEC_HEC_mod_change),row.names(AEC_HEC_no_mod_change))

##################################
#------get differential PSI events.
##################################
diff_path = "../Diff_PSI/diff_PSI_WC_test/pairwise_test/all_res/"
AEC_HEC_test_tab = read.table(file.path(diff_path, "AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_diff = AEC_HEC_test_tab[AEC_HEC_test_tab$FDR<0.05,]
AEC_HEC_not_diff = AEC_HEC_test_tab[AEC_HEC_test_tab$FDR>0.05,]

diff_type = c(rep("is_diff",nrow(AEC_HEC_diff)),rep("not_diff",nrow(AEC_HEC_not_diff)))
names(diff_type) = c(row.names(AEC_HEC_diff),row.names(AEC_HEC_not_diff))

all_event_type = data.frame(event_name = names(diff_type), diff_type)
all_event_type$event_change = event_change[as.character(all_event_type$event_name)]
all_event_type$mod_change_type = mod_change_type[as.character(all_event_type$event_name)]
all_event_type[is.na(all_event_type)] <-""
all_event_change_stat = as.data.frame(table(all_event_type[,2:3]))
all_event_change_filter = all_event_change_stat[all_event_change_stat$Freq>0,]
df = all_event_change_filter
df = df[order(df[,1],df[2],df[3]),]
#write.table(df,"modality_change_pie.txt",sep="\t",quote=F,row.names=F)
df$ymax = c(4381,4381+178,4381+1205+178,4381+1205+178+3440)
df$ymin = c(0,4381,4381+178,4381+1205+178)


#the order of Freq value should be from circle inside to outside.
library(ggplot2)
p <- ggplot(df) + 
#geom_rect(aes(fill = mod_change_type, ymax = ymax, ymin = ymin , xmax = 4, xmin = 3)) +
geom_rect(aes(fill = event_change, ymax = ymax , ymin = ymin , xmax = 3, xmin = 2)) +
geom_rect(aes(fill = diff_type, ymax = ymax , ymin = ymin , xmax = 2, xmin = 0)) +
xlim(c(0, 4))+ 
theme(aspect.ratio = 1)+
coord_polar(theta = 'y')
 
p
 
write.table(all_event_type,"./modality_diff_PSI/AEC_HEC_modality_diff_PSI_comparison.txt",sep="\t",row.names=F,quote=F)

##################################
#------get differential PSI events.
##################################
stage = c("AEC","HEC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])
stage_psi = event_psi[,samples]
library(data.table)
stage_psi_melt = melt(t(stage_psi))
colnames(stage_psi_melt) = c("sample","event_name","PSI")
stage_psi_melt_not_NA = na.omit(stage_psi_melt)

stage_psi_cv = t(apply(stage_psi,1,function(x){tapply(x,sample2stage[names(x),"stage"],function(x){sd(x,na.rm=T)})}))
stage_psi_cv[is.na(stage_psi_cv)]<-0
stage_psi_cv_AEC_high = stage_psi_cv[stage_psi_cv[,"AEC"]>0.1,]
stage_psi_cv_HEC_high = stage_psi_cv[stage_psi_cv[,"HEC"]>0.1,]
stage_psi_cv_high = union(row.names(stage_psi_cv_AEC_high),row.names(stage_psi_cv_HEC_high))


diff_path = "../Diff_PSI/diff_PSI_WC_test/pairwise_test/all_res/"
AEC_HEC_test_tab = read.table(file.path(diff_path, "AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_diff = AEC_HEC_test_tab[AEC_HEC_test_tab$FDR<0.05,]
AEC_HEC_not_diff = AEC_HEC_test_tab[AEC_HEC_test_tab$FDR>0.5,]
AEC_HEC_not_diff = AEC_HEC_not_diff[abs(AEC_HEC_not_diff$diff_PSI)<0.1,]

plot_PSI = stage_psi_melt_not_NA[stage_psi_melt_not_NA$event_name %in% row.names(AEC_HEC_not_diff),]
plot_PSI$stage = factor(as.character(sample2stage[plot_PSI$sample,"stage"]),levels = stage)
plot_PSI$diff = AEC_HEC_not_diff[as.character(plot_PSI$event_name),"diff_PSI"]
plot_PSI_high = plot_PSI[plot_PSI$event_name %in% stage_psi_cv_high, ]

plot_frm = plot_PSI
library(ggplot2)
library(ggpubr)
p = ggscatterhist(data = plot_frm, x = "diff", y = "PSI", shape=16,color="stage",
size=1,alpha=0.8,palette=c("#42a5f5","#9e9d24"),margin.plot="density",
margin.params=list(fill="stage",color="black",size=1),
legend=c(1,1),ggtheme = theme_minimal())

p






