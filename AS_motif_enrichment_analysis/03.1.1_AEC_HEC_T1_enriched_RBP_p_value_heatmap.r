rm(list=ls())
gc()

####################################
#get motif enrichment result
####################################
all_rMAPS_id = read.table("./rMAPS_enrichment/96_motif_RBP/rMAPS_motif_2_RBP_ID.txt",sep="\t",header=T,check.names=F)
row.names(all_rMAPS_id) = all_rMAPS_id$motif_id
#all_rMAPS_id$mouse_id = motif2gene_id[row.names(all_rMAPS_id),"gene_id"]
#write.table(all_rMAPS_id[,c(7,8,5,6,2:4,1)],"./rMAPS_enrichment/rMAPS_motif_2_RBP_ID.txt",sep="\t",quote=F,row.names=F)

#------get differential expressed_result
diff_path = "../../DEAnalysis/diff_res/gene_level/pairwise_test/all_res"
AEC_HEC_all = read.table(file.path(diff_path,"AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_diff = AEC_HEC_all[AEC_HEC_all$FDR<0.05,]
HEC_T1_all = read.table(file.path(diff_path,"HEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_diff = HEC_T1_all[HEC_T1_all$FDR<0.05,]
AEC_T1_all = read.table(file.path(diff_path,"AEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_diff = AEC_T1_all[AEC_T1_all$FDR<005,]

####################################
#AEC HEC p-value
####################################

AEC_HEC_enrich = read.table("./rMAPS_enrichment/HEC_T1_sig/HEC_T1_enriched_gene.txt",sep="\t",header=T)
row.names(AEC_HEC_enrich) = AEC_HEC_enrich[,1]
AEC_HEC_p = read.table("./rMAPS_enrichment/HEC_T1_sig/resultMaps/pVal.up.vs.bg.RNAmap.txt",sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_p_filter = AEC_HEC_p[row.names(AEC_HEC_enrich),c(3:6)]
AEC_HEC_log_p = -log10(AEC_HEC_p_filter)
colnames(AEC_HEC_log_p) = 1:4
AEC_HEC_p_tab = all_rMAPS_id[row.names(AEC_HEC_enrich),]
AEC_HEC_p_frm = data.frame(AEC_HEC_p_tab,AEC_HEC_diff[as.character(AEC_HEC_p_tab$mouse_id),])

library(data.table)
AEC_HEC_log_p_melt = melt(t(AEC_HEC_log_p))
colnames(AEC_HEC_log_p_melt) = c("region","motif_id","logP")
dot_plot_frm = AEC_HEC_log_p_melt
dot_plot_frm = na.omit(dot_plot_frm)
motif_sorted = c("SRp20.[AT]C[AT][AT]C","SRSF2.GGAG[AT][AGT]","SFPQ.[GT]T[AG][AG]T[GT][GT]",
"SRSF9.[GT]G[AG][AT]G[GC][AC]","SRSF9.A[GT]GA[ACG][AC][AG]","SRSF10.A[AG]AG[AG][AG][AG]")
dot_plot_frm$motif_id = factor(dot_plot_frm$motif_id,levels = rev(motif_sorted))
library(ggplot2)
library(RColorBrewer)
color =  brewer.pal(11,"RdYlBu")
color_Rd =  colorRampPalette(color[1:4])(100)
color_Yl =  colorRampPalette(color[4:8])(10)
color_Bu =  colorRampPalette(color[8:11])(30)
color = rev(c(color_Rd, color_Yl, color_Bu))
text_size = 8
line_size = 0.5/1.07

mytheme = theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA,color=NA),
strip.background = element_rect(fill=NA,color=NA),
strip.text = element_text(size=text_size,color="black"),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
panel.grid.minor = element_line(size=NA, color = NA),
axis.ticks = element_line(size=line_size,color="black"),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, hjust = 0.5, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')


p = ggplot(data = dot_plot_frm,aes(x = motif_id,y = region, color = logP))+
geom_point(shape=16,size=10)+
scale_color_gradientn(colors = color)+
mytheme+
coord_flip()
p


#------heatmap
library(pheatmap)
p = pheatmap(AEC_HEC_log_p,scale="none",cluster_col=F,cluster_row=F,color = color)


####################################
#HEC & T1 p-value
####################################
#------T1 upregulated
HEC_T1_enrich = read.table("./rMAPS_enrichment/HEC_T1_sig/HEC_T1_enriched_gene.txt",sep="\t",header=T)
row.names(HEC_T1_enrich) = HEC_T1_enrich[,1]
HEC_T1_up_p = read.table("./rMAPS_enrichment/HEC_T1_sig/resultMaps/pVal.up.vs.bg.RNAmap.txt",sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_up_p_filter = HEC_T1_up_p[row.names(HEC_T1_enrich),c(3:6)]
HEC_T1_up_log_p = -log10(HEC_T1_up_p_filter)
colnames(HEC_T1_up_log_p) = 1:4
HEC_T1_up_p_tab = all_rMAPS_id[row.names(HEC_T1_enrich),]
HEC_T1_up_p_frm = data.frame(HEC_T1_up_p_tab,HEC_T1_all[as.character(HEC_T1_up_p_tab$mouse_id),])

library(data.table)
HEC_T1_log_up_p_melt = melt(t(HEC_T1_up_log_p))
colnames(HEC_T1_log_up_p_melt) = c("region","motif_id","logP")
dot_plot_up_frm = HEC_T1_log_up_p_melt
dot_plot_up_frm$median_TPM = HEC_T1_up_p_frm[as.character(dot_plot_up_frm$motif_id),"T1_pre_HSC_median_TPM"]
dot_plot_up_frm = dot_plot_up_frm[dot_plot_up_frm$median_TPM>0,]
dot_plot_up_frm$type = rep("upregulated",nrow(dot_plot_up_frm))
up_order_tab = dot_plot_up_frm[order(dot_plot_up_frm$median_TPM,dot_plot_up_frm$motif_id),]

#------T1 downregulated
HEC_T1_down_p = read.table("./rMAPS_enrichment/AEC_HEC_T1/HEC_T1/resultMaps/pVal.dn.vs.bg.RNAmap.txt",sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_down_p_filter = HEC_T1_down_p[row.names(HEC_T1_enrich),c(3:6)]
HEC_T1_down_log_p = log10(HEC_T1_down_p_filter)
colnames(HEC_T1_down_log_p) = 1:4
HEC_T1_down_p_tab = all_rMAPS_id[row.names(HEC_T1_enrich),]
HEC_T1_down_p_frm = data.frame(HEC_T1_down_p_tab,HEC_T1_all[as.character(HEC_T1_down_p_tab$mouse_id),])

library(data.table)
HEC_T1_log_down_p_melt = melt(t(HEC_T1_down_log_p))
colnames(HEC_T1_log_down_p_melt) = c("region","motif_id","logP")
dot_plot_down_frm = HEC_T1_log_down_p_melt
dot_plot_down_frm$median_TPM = HEC_T1_down_p_frm[as.character(dot_plot_down_frm$motif_id),"HEC_median_TPM"]
dot_plot_down_frm = dot_plot_down_frm[dot_plot_down_frm$median_TPM>0,]
dot_plot_down_frm$type = rep("downregulated",nrow(dot_plot_down_frm))
down_order_tab = dot_plot_down_frm[order(dot_plot_down_frm$median_TPM,dot_plot_down_frm$motif_id),]
down_order_motif = unique(as.character(down_order_tab$motif_id))
order_samp = paste(rep(down_order_motif,each=2),c("upregulated","downregulated"),sep="_")

plot_frm = dot_plot_up_frm

values = seq(0,-log10(0.05)*5,by = -log10(0.05))
plot_frm_new = plot_frm$logP
names(plot_frm_new)
plot_frm_new[which(abs(plot_frm_new)<(-log10(0.05)))]<-NA
plot_frm$new_logP = plot_frm_new

library(ggplot2)
library(RColorBrewer)
color =  brewer.pal(11,"RdYlBu")
color_Rd =  colorRampPalette(color[1:3])(100)
color_Yl =  colorRampPalette(color[3:9])(50)
color_Bu =  colorRampPalette(color[9:11])(100)

#colors = rev(colorRampPalette(color)(100))
colors = rev(c(color_Rd, color_Yl, color_Bu))
text_size = 8
line_size = 0.5/1.07
p = ggplot(data = plot_frm,aes(x = motif_id,y = region, color = new_logP,size = median_TPM))+
geom_point(shape=16)+
scale_color_gradientn(colors = colors,limits = c(log10(0.05)*5,-log10(0.05)*5),breaks=seq(-6,6,by=3))+
coord_flip()
p


####################################
#T1 only included p-value
####################################
#------T1 upregulated
HEC_T1_enrich = read.table("./rMAPS_enrichment/T1_only_sig/T1_only_enriched_gene.txt",sep="\t",header=T)
row.names(HEC_T1_enrich) = HEC_T1_enrich[,1]
HEC_T1_p = read.table("./rMAPS_enrichment/T1_only_sig/resultMaps/pVal.up.vs.bg.RNAmap.txt",sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_p_filter = HEC_T1_p[row.names(HEC_T1_enrich),c(3:6)]
HEC_T1_log_p = -log10(HEC_T1_p_filter)
colnames(HEC_T1_log_p) = 1:4
HEC_T1_p_tab = all_rMAPS_id[row.names(HEC_T1_enrich),]
HEC_T1_p_frm = data.frame(HEC_T1_p_tab,HEC_T1_all[as.character(HEC_T1_p_tab$mouse_id),])
HEC_T1_median_TPM = apply(HEC_T1_p_frm,1,function(x){median(as.numeric(x[c("HEC_median_TPM","T1_pre_HSC_median_TPM")]))})

library(data.table)
HEC_T1_log_p_melt = melt(t(HEC_T1_log_p))
colnames(HEC_T1_log_p_melt) = c("region","motif_id","logP")
dot_plot_frm = HEC_T1_log_p_melt
dot_plot_frm$median_TPM = HEC_T1_median_TPM[as.character(dot_plot_frm$motif_id)]
dot_plot_frm = dot_plot_frm[dot_plot_frm$median_TPM>0,]
dot_plot_frm = dot_plot_frm[order(dot_plot_frm$median_TPM,dot_plot_frm$motif_id),]
sort_samp = unique(as.character(dot_plot_frm$motif_id))

plot_frm = dot_plot_frm
plot_frm$motif_id= factor(plot_frm$motif_id,levels=sort_samp)
summary(plot_frm$logP)
values = seq(0,-log10(0.05)*4,by = -log10(0.05))
plot_frm_new = plot_frm$logP
plot_frm_new[which(abs(plot_frm_new)<(-log10(0.05)))]<-NA
plot_frm$new_logP = plot_frm_new

library(ggplot2)
library(RColorBrewer)
color =  brewer.pal(11,"RdYlBu")
color_Rd =  colorRampPalette(color[1:5])(100)
color_Yl =  colorRampPalette(color[3:9])(50)
color_Bu =  colorRampPalette(color[9:11])(100)

colors = rev(color_Rd)
#colors = rev(c(color_Rd, color_Yl, color_Bu))
text_size = 8
line_size = 0.5/1.07
p = ggplot(data = plot_frm,aes(x = motif_id,y = region, color = new_logP,size = median_TPM))+
geom_point(shape=16)+
scale_color_gradientn(colors = colors,limits = c(0,4.2),breaks=seq(0,4,by=2))+
coord_flip()
p









