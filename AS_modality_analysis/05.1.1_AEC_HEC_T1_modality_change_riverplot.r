rm(list=ls())
gc()
#install.packages("riverplot")
library(riverplot)
library(data.table)
###########################################################
#load data
###########################################################
events_modality = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1)

get_mod_tab <- function(stage1,stage2, mod_tab = events_modality){
stage_modality = mod_tab[,c(stage1,stage2)]
stage_modality_filter = stage_modality[apply(stage_modality,1,function(x){length(x[x==""])<length(c(stage1,stage2))}),]
stage_modality_stat = table(stage_modality)
row.names(stage_modality_stat)[row.names(stage_modality_stat)==""]="no_modality"
colnames(stage_modality_stat)[colnames(stage_modality_stat)==""]="no_modality"
sort_name = c("no_modality","multimodal", "bimodal", "middle", "excluded", "included")
stage_modality_stat_sorted = stage_modality_stat[intersect(sort_name,row.names(stage_modality_stat)),
intersect(sort_name,colnames(stage_modality_stat))]
row.names(stage_modality_stat_sorted) = paste(stage1,row.names(stage_modality_stat_sorted),sep="_")
colnames(stage_modality_stat_sorted) = paste(stage2,colnames(stage_modality_stat_sorted),sep="_")
stage_modality_melt = melt(stage_modality_stat_sorted)
colnames(stage_modality_melt) = c("N1","N2","Value")
return(stage_modality_melt)
}

AEC_HEC_stat = get_mod_tab(stage1 = "AEC", stage2 = "HEC")
HEC_T1_stat = get_mod_tab(stage1 = "HEC", stage2 = "T1_pre_HSC")
T1_T2_stat = get_mod_tab(stage1 = "T1_pre_HSC", stage2 = "T2_pre_HSC")
T2_E12_stat = get_mod_tab(stage1 = "T2_pre_HSC", stage2 = "E12")
E12_E14_stat = get_mod_tab(stage1 = "E12", stage2 = "E14")
E14_Adult_stat = get_mod_tab(stage1 = "E14", stage2 = "Adult_HSC")

AEC_HEC_T1_frm = rbind(AEC_HEC_stat, HEC_T1_stat)
###########################################################
#riverplot
###########################################################

data_edges = AEC_HEC_T1_frm
edges_filter = data_edges[data_edges$Value>0,]
edges_filter$Value = log2(edges_filter$Value+1)
data_nodes = data.frame(ID = unique(c(as.character(data_edges$N1), as.character(data_edges$N2))), stringsAsFactors = FALSE)  
data_nodes$x = rep(c(1:3),each=6)
data_nodes$y = rep(c(0:5),n=3)
rownames(data_nodes) = data_nodes$ID


library(RColorBrewer)

palette = c("#808080","#E79F58","#9F58E7","#9FE758","#589FE7","#E75858")

styles = lapply(data_nodes$y, function(n) {
  list(col = palette[n+1], lty = 0, textcol = "black")
})
names(styles) = data_nodes$ID

rp <- list(nodes = data_nodes, edges = edges_filter, styles = styles)
class(rp) <- c(class(rp), "riverplot")
plot(rp, plot_area = 1, yscale=0.02)






