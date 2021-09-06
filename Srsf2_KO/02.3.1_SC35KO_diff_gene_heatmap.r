rm(list=ls())
gc()
indir = "../expression_data/download_from_linux/salmon_merge_counts"
gene_tpm = read.table(paste(indir,"gene_normalized_TPM.txt",sep="/"),sep="\t",header=T,row.names=1,check.names=F)

diff_res = read.table("./diff_res/EC_WT_KO_edgeR_test_sig_res.txt",sep="\t",header=T,row.names=1)

diff_gene_exp = gene_tpm[intersect(row.names(gene_tpm),row.names(diff_res)),]
diff_gene_exp_filter = diff_gene_exp[rowMeans(diff_gene_exp)>1,]
plot_frm = log2(diff_gene_exp_filter+1)

library(RColorBrewer)
color =  brewer.pal(11,"RdYlBu")
color_Rd =  colorRampPalette(color[1:4])(100)
color_Yl =  colorRampPalette(color[4:8])(50)
color_Bu =  colorRampPalette(color[8:11])(100)
color = rev(c(color_Rd, color_Yl, color_Bu))

library(pheatmap)
q = pheatmap(plot_frm,scale="row",show_rownames=F,cutree_row=2,cutree_col=2)

candidates_gene = read.table("./diff_res/SC35_KO_diff_gene_candidates.txt",sep="\t",header=T,row.names=1)
candidates_gene_name = as.character(candidates_gene$gene_name)
candidates_row = data.frame(gene_name = candidates_gene_name)
row.names(candidates_row) = row.names(candidates_gene)
other_gene = data.frame(gene_name=rep('',length(setdiff(row.names(plot_frm),row.names(candidates_gene)))))
row.names(other_gene) = setdiff(row.names(plot_frm),row.names(candidates_gene))
anno_row = rbind(other_gene, candidates_row)

p = pheatmap(plot_frm,scale="row",show_rownames=T,cutree_row=2,cutree_col=2, labels_row = anno_row[row.names(plot_frm),1])
