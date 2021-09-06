rm(list=ls())
gc()

sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T,stringsAsFactors=F)
row.names(sample2stage) = sample2stage$sample
stage2color = read.table("../../sample_infor/stage2color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)


event_psi = read.table("../merge_PSI/merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,row.names=1,check.names=F)
modality_tab = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1)

AEC_HEC_in = read.table("./AEC_HEC_T1/AEC_HEC_included_events.txt",sep="\t",header=T,check.names=F)
HEC_T1_in = read.table("./AEC_HEC_T1/HEC_T1_included_events.txt",sep="\t",header=T,check.names=F)
AEC_T1_in = read.table("./AEC_HEC_T1/AEC_T1_included_events.txt",sep="\t",header=T,check.names=F)
AEC_HEC_T1_in_events = Reduce(union,list(row.names(AEC_HEC_in),row.names(HEC_T1_in),row.names(AEC_T1_in)))

stages = c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stages,])
T1_sig_psi = event_psi[AEC_HEC_T1_in_events,samples]

library(RColorBrewer)
color =  brewer.pal(11,"RdYlBu")
color_Rd =  colorRampPalette(color[1:4])(100)
color_Yl =  colorRampPalette(color[4:8])(50)
color_Bu =  colorRampPalette(color[8:11])(100)
color = rev(c(color_Rd, color_Yl, color_Bu))

anno_col = data.frame(type = sample2stage[samples,'stage'])
row.names(anno_col) = samples
anno_col$type=factor(anno_col$type,levels=stages)


plot_frm = T1_sig_psi
plot_frm[is.na(plot_frm)]=0
library(pheatmap)
cluster_num = 4
p = pheatmap(plot_frm,scale="row",col = color,cluster_col=F,cluster_row=T,show_rownames=F,
cutree_row=cluster_num,cutree_col=1,show_colnames=F,clustering_distance_rows = "correlation",
gaps_col = c(18,42))
row_order = p$tree_row$order
#------how to change cluster order
hclust.result<-p$tree_row
cluster = cutree(hclust.result,k=cluster_num)
cluster_name = paste("C",cluster,sep="")
names(cluster_name) = names(cluster)
hclust.result.order<-p$tree_row$order
cluster_name_ordered = cluster_name[hclust.result.order]
cluster_categroy = unique(cluster_name_ordered)
new_cluster_name = paste("cluster",1:length(cluster_categroy),sep="")
names(new_cluster_name) = cluster_categroy

row_cluster_name = data.frame(old_cluster = cluster_name)
row_cluster_name$new_cluster = new_cluster_name[as.character(row_cluster_name$old_cluster)]
cluster_frm = T1_sig_psi[row_order,]
cluster_frm$cluster = row_cluster_name[row.names(cluster_frm),"new_cluster"]
cluster_tab = data.frame(cluster_frm,modality_tab[row.names(cluster_frm),stages])
write.table(cluster_tab,"./figure/AEC_HEC_T1_included_events_heatmap.txt",sep="\t",row.names=F,quote=F)


anno_row = data.frame(cluster = row_cluster_name[,2])
row.names(anno_row) = row.names(row_cluster_name)
plot_new = T1_sig_psi[row.names(cluster_tab),]
p = pheatmap(plot_new,scale="none",cluster_col=F,cluster_row=F,show_rownames=F,
show_colnames=F,gaps_col = c(18,42),annotation_col = anno_col,annotation_row = anno_row)

col = color,

#------each cluster violin plot.
cluster_PSI_tab = data.frame(cluster_tab[,c("event_name","cluster")],T1_sig_psi[row.names(cluster_tab),],check.names=F)
cluster_PSI_median = apply(cluster_PSI_tab[,-c(1:2)],2,function(x){tapply(x,as.character(cluster_PSI_tab$cluster),function(y){median(y,na.rm=T)})})
library(data.table)
cluster_PSI_melt = melt(cluster_PSI_median)
colnames(cluster_PSI_melt) = c("cluster","sample","PSI")
cluster_PSI_melt$stage = as.character(sample2stage[as.character(cluster_PSI_melt$sample),"stage"])
cluster_PSI_melt$stage = factor(cluster_PSI_melt$stage,levels = stages)

eight_stage_col = stage2color[stages,"figurecolor"]
plot_frm = cluster_PSI_melt
text_size = 8
line_size = 0.5/1.07
library(ggplot2)
p <- ggplot(plot_frm,mapping = aes(x = stage, y = PSI,color = stage))+
geom_violin(position = "dodge", na.rm = TRUE, fill=NA, scale = "width", width = 0.8)+
#geom_boxplot(fill=NA,width=0.3,outlier.shape=NA)+
geom_jitter(aes(x = stage, y = PSI, color = PSI),width=0.1,size=0.5,shape=16,color="black")+
scale_color_manual(values = eight_stage_col)+
scale_y_continuous(limits = c(-0.01,1.01),breaks = seq(0,1,by=0.5))+
labs(x="",y='PSI')+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_wrap(~cluster,ncol=1)

p

library(plyr)
line_plot_frm = plot_frm
line_plot_frm = ddply(line_plot_frm,.(cluster,stage),transform,median_PSI=median(PSI,na.rm=T))
line_plot_frm_uniq = unique(line_plot_frm[,c(1,4:5)])
line_plot_frm_uniq$index = rep(1:8,5)
q <- ggplot(line_plot_frm_uniq,mapping = aes(x = index, y = median_PSI))+
#geom_violin(position = "dodge", na.rm = TRUE, fill=NA, scale = "width", width = 0.8)+
#geom_boxplot(fill=NA,width=0.3,outlier.shape=NA)+
#geom_jitter(width=0.1,size=0.5,shape=16,color="black")+
geom_point()+
stat_smooth(method = 'loess', color="black",size=line_size*1.5, n= 8, span = 0.3, se = F)+
scale_color_manual(values = eight_stage_col)+
scale_y_continuous(limits = c(-0.01,1.01),breaks = seq(0,1,by=0.5))+
labs(x="",y='PSI')+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_wrap(~cluster,ncol=1)

q










cluster_tab_modality = event_modality[row.names(cluster_tab),stages]
cluster_tab_mod_stat = apply(cluster_tab_modality,2,function(x){
tapply(x,row_cluster_name[names(x),"new_cluster"],function(y){length(x[which(x=="included")])})
})
all_num = table(cluster_tab$cluster)
EC_YS_in = table(cluster_tab$cluster,cluster_tab$EC_YS)[,"included"]
EMP_in = table(cluster_tab$cluster,cluster_tab$EMP)[,"included"]
EC_in = table(cluster_tab$cluster,cluster_tab$EC)[,"included"]
T1_in = table(cluster_tab$cluster,cluster_tab$T1_pre_HSC)[,"included"]
T2_in = table(cluster_tab$cluster,cluster_tab$T2_pre_HSC)[,"included"]
E12_in = table(cluster_tab$cluster,cluster_tab$E12)[,"included"]
E14_in = table(cluster_tab$cluster,cluster_tab$E14)[,"included"]
Adult_in = table(cluster_tab$cluster,cluster_tab$Adult_HSC)[,"included"]
included_num = data.frame(all_num, EC_YS = EC_YS_in, EMP = EMP_in, EC = EC_in, T1_pre_HSC = T1_in,
T2_pre_HSC = T2_in, E12 = E12_in,E14 = E14_in,Adult_HSC = Adult_in)
included_num_tab = included_num[,-1]
included_num_propor = t(apply(included_num_tab,1,function(x){x[-1]/x[1]}))


pheatmap(included_num_propor,scale="none",cluster_col=F)
library(data.table)
propor_melt = melt(included_num_propor)
colnames(propor_melt) = c("cluster","stage","propor")

ggplot_color = colorRampPalette(c("#4292C6","#EF3B2C"))(100)
library(ggplot2)
p = ggplot(propor_melt, aes(x = stage,y=propor,fill = propor))+
geom_bar(stat="identity",width=.7)+
scale_fill_gradientn(colors = ggplot_color)+
facet_wrap(.~cluster,ncol=1)
p




plot_frm_order = T1_sig_psi[row_order,]
plot_frm_scale = t(apply(plot_frm_order,1,scale))
colnames(plot_frm_scale) = colnames(plot_frm_order)

anno_row = data.frame(cluster = row_cluster_name$new_cluster)
row.names(anno_row) = row.names(row_cluster_name)
p = pheatmap(plot_frm_order,scale="none",cluster_col=F,cluster_row=F,show_rownames=F,
show_colnames=F,gaps_col = c(27,53,70,96,134,155,187),gaps_row = c(85,308,386,585),
annotation_row = anno_row,annotation_col=anno_col)

















