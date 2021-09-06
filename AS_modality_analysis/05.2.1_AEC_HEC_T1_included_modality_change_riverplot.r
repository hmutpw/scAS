rm(list=ls())
gc()
#install.packages("ggforce")
library(ggforce)
library(ggplot2)
library(data.table)
###########################################################
#load data
###########################################################
stage2color = read.table("../../sample_infor/all_color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)
events_modality = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1)

###########################################################
#------all events riverplot
###########################################################
stage_all = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
stage_modality = events_modality[,stage_all] 
stage_modality_filter = stage_modality[apply(stage_modality,1,function(x){length(x[which(x!="")])>0}),]
sort_name = c("no_modality","multimodal", "bimodal", "middle", "excluded", "included")
all_events_mat = as.matrix(stage_modality_filter)
all_events_mat[all_events_mat==""] <- "no_modality"
all_events_mod_stat = as.data.frame(table(as.data.frame(all_events_mat)))

data_edges = all_events_mod_stat
edges_filter = data_edges[data_edges$Freq>0,]

data = edges_filter
data <- gather_set_data(data, 1:7)
data$x = factor(data$x,levels=stage_all)
data$y = factor(data$y,levels=rev(sort_name))
data$T1_pre_HSC = factor(data$T1_pre_HSC,levels=rev(sort_name))

stage_color = stage2color[c("red","cyanblue","green","purple","orange","gray"),"num4"]

p = ggplot(data, aes(x, id = id, split = y, value = Freq)) +
geom_parallel_sets(aes(fill = T1_pre_HSC), alpha = 0.5, axis.width = 0.1) +
geom_parallel_sets_axes(axis.width = 0.1) +
geom_parallel_sets_labels(colour = 'white')+
scale_color_manual(values=stage_color)+
scale_fill_manual(values=stage_color)

p
###########################################################
#------AEC HEC T1 events riverplot
###########################################################
stage = c("AEC","HEC","T1_pre_HSC")

AEC_HEC_T1_mod_tab = events_modality[,stage]
AEC_HEC_T1_included = AEC_HEC_T1_mod_tab[apply(AEC_HEC_T1_mod_tab,1,function(x){length(x[which(x!="")])>0}),]
AEC_HEC_T1_included_mat = as.matrix(AEC_HEC_T1_included)
AEC_HEC_T1_included_mat[AEC_HEC_T1_included_mat==""] <- "no_modality"
AEC_HEC_T1_included_stat = as.data.frame(table(as.data.frame(AEC_HEC_T1_included_mat)))

data_edges = AEC_HEC_T1_included_stat
edges_filter = data_edges[data_edges$Freq>0,]


data = edges_filter
data <- gather_set_data(data, 1:length(stage))
data$x = factor(data$x,levels=stage)
data$y = factor(data$y,levels=rev(sort_name))
data$T1_pre_HSC = factor(data$T1_pre_HSC,levels=rev(sort_name))

stage_color = stage2color[c("red","cyanblue","green","purple","orange","gray"),"num4"]

p = ggplot(data, aes(x, id = id, split = y, value = Freq)) +
geom_parallel_sets(fill = "red", alpha = 0.5, axis.width = 0.1) +
geom_parallel_sets_axes(axis.width = 0.1) +
geom_parallel_sets_labels(colour = 'white')+
scale_color_manual(values=stage_color)+
scale_fill_manual(values=stage_color)

p


###########################################################
#------T1 included events riverplot
###########################################################
stage = c("AEC","HEC","T1_pre_HSC")

AEC_HEC_T1_mod_tab = events_modality[,stage]
AEC_HEC_T1_included = AEC_HEC_T1_mod_tab[AEC_HEC_T1_mod_tab$T1_pre_HSC=="included",]
AEC_HEC_T1_included_mat = as.matrix(AEC_HEC_T1_included)
AEC_HEC_T1_included_mat[AEC_HEC_T1_included_mat==""] <- "no_modality"
AEC_HEC_T1_included_stat = as.data.frame(table(as.data.frame(AEC_HEC_T1_included_mat)))

data_edges = AEC_HEC_T1_included_stat
edges_filter = data_edges[data_edges$Freq>0,]


data = edges_filter
data <- gather_set_data(data, 1:length(stage))
data$x = factor(data$x,levels=stage)
data$y = factor(data$y,levels=rev(sort_name))
data$HEC = factor(data$HEC,levels=rev(sort_name))

stage_color = stage2color[c("red","cyanblue","green","orange","gray"),"num4"]

p = ggplot(data, aes(x, id = id, split = y, value = Freq)) +
geom_parallel_sets(aes(fill = HEC), alpha = 0.5, axis.width = 0.1) +
geom_parallel_sets_axes(axis.width = 0.1) +
geom_parallel_sets_labels(colour = 'white')+
scale_color_manual(values=stage_color)+
scale_fill_manual(values=stage_color)

p



