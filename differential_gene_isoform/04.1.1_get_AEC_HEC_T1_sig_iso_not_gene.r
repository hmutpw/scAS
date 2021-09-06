rm(list=ls())
gc()

get_signature <- function(pos_stage, neg_stage, exp_cutoff = 1, in_dir){
stagepare = paste(neg_stage, pos_stage, sep="_")
#------isoform level.
iso_files = file.path(in_dir, "transcript_level/pairwise_test/all_res",paste(stagepare,"_test_all_res.txt",sep=""))
iso_tab = read.table(iso_files,sep="\t",header=T,row.names=1,stringsAsFactors=F)
iso_high = iso_tab[apply(iso_tab[,grep("_median_TPM$",colnames(iso_tab))],1,function(x){length(x[which(x>exp_cutoff)])>0}),]
diff_iso = iso_high[iso_high$FDR<0.05,]
diff_iso = na.omit(diff_iso)
diff_iso_up = diff_iso[diff_iso$logFC>0,]
diff_iso_down = diff_iso[diff_iso$logFC<0,]

#------gene level.
gene_files = file.path(in_dir, "gene_level/pairwise_test/all_res",paste(stagepare,"_test_all_res.txt",sep=""))
gene_tab = read.table(gene_files,sep="\t",header=T,row.names=1,stringsAsFactors=F)
gene_high = gene_tab[apply(gene_tab[,grep("_median_TPM$",colnames(gene_tab))],1,function(x){length(x[which(x>exp_cutoff)])>0}),]
diff_gene = gene_high[gene_high$FDR<0.05,]
diff_gene = na.omit(diff_gene)
diff_gene_up = diff_gene[diff_gene$logFC>0,]
diff_gene_down = diff_gene[diff_gene$logFC<0,]
not_diff_gene = gene_tab[gene_tab$FDR>0.1,]

gene_tab_out = gene_tab[,-c(1:2)]
colnames(gene_tab_out) = paste("gene_",colnames(gene_tab_out),sep="")
#------up iso not gene.
pos_up_iso_gene_id = unique(as.character(diff_iso_up$gene_id))
pos_up_iso_not_gene = intersect(pos_up_iso_gene_id,union(row.names(not_diff_gene),row.names(diff_gene_down)))
pos_up_iso_not_gene_iso_id = row.names(diff_iso_up[diff_iso_up$gene_id %in% pos_up_iso_not_gene,])
pos_up_iso_not_gene_iso_type = rep(paste(pos_stage,"_enrich",sep=""),length(pos_up_iso_not_gene_iso_id))
names(pos_up_iso_not_gene_iso_type) = pos_up_iso_not_gene_iso_id

#------down iso not gene.
pos_down_iso_gene_id = unique(as.character(diff_iso_down$gene_id))
pos_down_iso_not_gene = intersect(pos_down_iso_gene_id,union(row.names(not_diff_gene),row.names(diff_gene_up)))
pos_down_iso_not_gene_iso_id = row.names(diff_iso_down[diff_iso_down$gene_id %in% pos_down_iso_not_gene,])
pos_down_iso_not_gene_iso_type = rep(paste(neg_stage,"_enrich",sep=""),length(pos_down_iso_not_gene_iso_id))
names(pos_down_iso_not_gene_iso_type) = pos_down_iso_not_gene_iso_id

#------merge two types of isoforms.
pos_neg_iso_not_gene_iso_type = c(pos_up_iso_not_gene_iso_type, pos_down_iso_not_gene_iso_type)
pos_neg_iso_not_gene_iso_id = names(pos_neg_iso_not_gene_iso_type)
pos_neg_iso_not_gene_iso_pair = rep(paste(neg_stage, pos_stage, sep="_"),length(pos_neg_iso_not_gene_iso_id))
names(pos_neg_iso_not_gene_iso_pair) = pos_neg_iso_not_gene_iso_id

#------get sig iso not gene all iso table.
pos_neg_iso_not_gene_iso_gene_id = unique(as.character(iso_high[pos_neg_iso_not_gene_iso_id,"gene_id"]))
pos_neg_iso_not_gene_iso_all_tab = iso_high[iso_high$gene_id %in% pos_neg_iso_not_gene_iso_gene_id,]
pos_neg_iso_not_gene_iso_all_out = data.frame(transcript_id = row.names(pos_neg_iso_not_gene_iso_all_tab),
pos_neg_iso_not_gene_iso_all_tab, gene_tab_out[as.character(pos_neg_iso_not_gene_iso_all_tab$gene_id),],
enrich_type = pos_neg_iso_not_gene_iso_type[row.names(pos_neg_iso_not_gene_iso_all_tab)],
stagepare = pos_neg_iso_not_gene_iso_pair[row.names(pos_neg_iso_not_gene_iso_all_tab)],check.names=F)
out_frm = pos_neg_iso_not_gene_iso_all_out
out_frm_sorted = out_frm[order(out_frm$gene_name,out_frm$FDR),]

return(out_frm_sorted)
}

out_path = "./diff_iso_not_gene/AEC_HEC_T1"
dir.create(out_path,showWarnings = FALSE, recursive = T)
AEC_HEC_sig = get_signature(pos_stage = "HEC", neg_stage = "AEC", in_dir="./diff_res")
write.table(AEC_HEC_sig,file.path(out_path,"AEC_HEC_diff_iso_not_gene.txt"),sep="\t",row.names=F,quote=F)

HEC_T1_sig = get_signature(pos_stage = "T1_pre_HSC", neg_stage = "HEC", in_dir="./diff_res")
HEC_T1_T1_enrich_gene = unique(as.character(HEC_T1_sig[HEC_T1_sig$enrich_type=="T1_pre_HSC_enrich","gene_id"]))
HEC_T1_sig_new = HEC_T1_sig[HEC_T1_sig$gene_id %in% HEC_T1_T1_enrich_gene,]
write.table(HEC_T1_sig_new,file.path(out_path,"HEC_T1_pre_HSC_diff_iso_not_gene.txt"),sep="\t",row.names=F,quote=F)

AEC_T1_sig = get_signature(pos_stage = "T1_pre_HSC", neg_stage = "AEC", in_dir="./diff_res")
AEC_T1_T1_enrich_gene = unique(as.character(AEC_T1_sig[AEC_T1_sig$enrich_type=="T1_pre_HSC_enrich","gene_id"]))
AEC_T1_sig_new = AEC_T1_sig[AEC_T1_sig$gene_id %in% AEC_T1_T1_enrich_gene,]
write.table(AEC_T1_sig_new,file.path(out_path,"AEC_T1_pre_HSC_diff_iso_not_gene.txt"),sep="\t",row.names=F,quote=F)

AEC_HEC_T1_diff_iso_not_gene = rbind(AEC_HEC_sig[,c(1:6,9:11,14:18)], HEC_T1_sig_new[,c(1:6,9:11,14:18)])
write.table(AEC_HEC_T1_diff_iso_not_gene,file.path(out_path,"AEC_HEC_T1_pre_HSC_diff_iso_not_gene_with_all_iso.txt"),sep="\t",row.names=F,quote=F)

AEC_HEC_T1_diff_iso_not_gene_out = AEC_HEC_T1_diff_iso_not_gene[!is.na(AEC_HEC_T1_diff_iso_not_gene$enrich_type),]
write.table(AEC_HEC_T1_diff_iso_not_gene_out,file.path(out_path,"AEC_HEC_T1_pre_HSC_diff_iso_not_gene_merge.txt"),sep="\t",row.names=F,quote=F)

out_path = "./diff_iso_not_gene/T1_T2"
dir.create(out_path,showWarnings = FALSE, recursive = T)
T1_T2_sig = get_signature(pos_stage = "T2_pre_HSC", neg_stage = "T1_pre_HSC", in_dir="./diff_res")
write.table(T1_T2_sig,file.path(out_path,"T1_T2_diff_iso_not_gene.txt"),sep="\t",row.names=F,quote=F)








