rm(list=ls())
gc()
gene_tpm = read.table("../../expression_data/merged_exp/7_stage_gene_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage=c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])

#------get differential test result.
test_path = "../../DEAnalysis/DESeq2_res/gene_level/all_res"
AEC_HEC_res = read.table(file.path(test_path,"AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_res = read.table(file.path(test_path,"HEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_res = read.table(file.path(test_path,"AEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)

diff_path = "../../DEAnalysis/DESeq2_res/gene_level/diff_res"
AEC_HEC_diff = read.table(file.path(diff_path,"AEC_HEC_test_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_diff = read.table(file.path(diff_path,"HEC_T1_pre_HSC_test_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_diff = read.table(file.path(diff_path,"AEC_T1_pre_HSC_test_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_T1_diff_gene = Reduce(union,list(row.names(AEC_HEC_diff), row.names(HEC_T1_diff),row.names(AEC_T1_diff)))

#------get rMAPS RBP gene.
RBP_tab = read.table("./rMAPS_enrichment/96_motif_RBP/rMAPS_motif_2_RBP_ID.txt",sep="\t",header=T,row.names=1,check.names=F)
RBP_gene_id = unique(as.character(RBP_tab$mouse_id))

AEC_HEC_T1_diff_RBP_gene = intersect(AEC_HEC_T1_diff_gene, RBP_gene_id)
AEC_HEC_T1_diff_RBP_gene_tpm = gene_tpm[AEC_HEC_T1_diff_RBP_gene,samples]

AEC_HEC_diff_RBP_gene = AEC_HEC_diff[intersect(row.names(AEC_HEC_diff),RBP_gene_id),]
HEC_T1_diff_RBP_gene = HEC_T1_diff[intersect(row.names(HEC_T1_diff),RBP_gene_id),]
AEC_T1_diff_RBP_gene = AEC_T1_diff[intersect(row.names(AEC_T1_diff),RBP_gene_id),]



#------plot heatmap.

plot_frm = log2(AEC_HEC_T1_diff_RBP_gene_tpm+1)

cluster_num = 1
library(pheatmap)
p = pheatmap(plot_frm,scale="row",cutree_row=cluster_num,cluster_col=T,
show_rownames=F,fontsize = 8)
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
anno_row = data.frame(cluster = cluster_name_ordered)
anno_row$cluster = new_cluster_name[as.character(anno_row$cluster)]
anno_row$cluster = factor(anno_row$cluster,
levels = paste("cluster",1:cluster_num,sep=""))

cluster_frm = data.frame(gene_id = row.names(anno_row),
gene_tpm[row.names(anno_row),c("gene_name","gene_type")],
cluster = anno_row[,"cluster"],
candidate_gene_median[row.names(anno_row),colnames(plot_frm)],check.names=F)
out_dir = "./WT_D2_D4_dynamic/heatmap_clusters"
dir.create(out_dir,recursive = T)
write.table(cluster_frm,file.path(out_dir,"WT_D2_D4_input_not_diff_loc_diff_gene_heatmap_clusters.txt"),sep="\t",quote=F,row.names=F)

plot_frm_ordered = plot_frm[row.names(cluster_frm),]
p = pheatmap(plot_frm_ordered,scale="row",cluster_row=F,cluster_col=F,
show_rownames=F,annotation_row = anno_row, treeheight_row=F,fontsize = 8,
gaps_col = seq(3,ncol(plot_frm),by=3), main="Input_not_diff_other_locs_diff",
legend_breaks=seq(-3,3,by=3))













