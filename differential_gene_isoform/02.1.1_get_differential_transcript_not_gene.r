rm(list=ls())
gc()

get_signature <- function(stagepare, exp_cutoff = 1, in_dir, out_dir){
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
#not_diff_gene = gene_tab[gene_tab$FDR>0.1,]

#------get sig iso & gene
sig_iso_id = row.names(diff_iso_up)
sig_iso_gene_id = unique(as.character(diff_iso_up$gene_id))
sig_gene = row.names(diff_gene_up)
sig_iso_not_gene = setdiff(sig_iso_gene_id, sig_gene)
sig_iso_gene = intersect(sig_iso_gene_id, sig_gene)

sig_iso_not_gene_iso_id = row.names(diff_iso_up[diff_iso_up$gene_id %in% sig_iso_not_gene,])
sig_iso_gene_iso_id = row.names(diff_iso_up[diff_iso_up$gene_id %in% sig_iso_gene,])
sig_type = c(rep("isoform_only",length(sig_iso_not_gene_iso_id)),
rep("gene_isoform",length(sig_iso_gene_iso_id)))
names(sig_type) = c(sig_iso_not_gene_iso_id,sig_iso_gene_iso_id)

#------add non sig isoform to output.
sig_iso_sig_gene_id = unique(as.character(iso_tab[sig_iso_gene_iso_id,"gene_id"]))
sig_iso_sig_gene_tab = iso_tab[iso_tab$gene_id %in% sig_iso_sig_gene_id,]
sig_iso_sig_gene_out_frm = data.frame(transcript_id = row.names(sig_iso_sig_gene_tab),
sig_iso_sig_gene_tab, gene_tab[as.character(sig_iso_sig_gene_tab$gene_id),-c(1:3)],
sig_type = sig_type[row.names(sig_iso_sig_gene_tab)], check.names=F)

sig_iso_gene_path = file.path(out_dir,"sig_iso_sig_gene")
dir.create(sig_iso_gene_path,recursive = T)
write.table(sig_iso_sig_gene_out_frm,file.path(sig_iso_gene_path,paste(stagepare,"_sig_iso_sig_gene_res.txt",sep="")),
sep="\t",row.names=F,quote=F)

#------add sig iso not sig gene
sig_iso_not_sig_gene_id = unique(as.character(iso_tab[sig_iso_not_gene_iso_id,"gene_id"]))
sig_iso_not_sig_gene_tab = iso_tab[iso_tab$gene_id %in% sig_iso_not_sig_gene_id,]
sig_iso_not_sig_gene_out_frm = data.frame(transcript_id = row.names(sig_iso_not_sig_gene_tab),
sig_iso_not_sig_gene_tab,gene_tab[as.character(sig_iso_not_sig_gene_tab$gene_id),-c(1:3)],
sig_type = sig_type[row.names(sig_iso_not_sig_gene_tab)], check.names=F)

sig_iso_not_gene_path = file.path(out_dir,"sig_iso_not_gene")
dir.create(sig_iso_not_gene_path,recursive = T)
write.table(sig_iso_not_sig_gene_out_frm,file.path(sig_iso_not_gene_path,paste(stagepare,"_sig_iso_not_gene_res.txt",sep="")),
sep="\t",row.names=F,quote=F)

sig_iso_frm = list(sig_iso_gene = sig_iso_sig_gene_out_frm,sig_iso_only = sig_iso_not_sig_gene_out_frm)
return(sig_iso_frm)
}
out_path = "./AEC_HEC_T1_sig"
AEC_HEC_sig = get_signature(stagepare="AEC_HEC", in_dir="./diff_res",out_dir = out_path)

HEC_T1_sig = get_signature(stagepare="HEC_T1_pre_HSC", in_dir="./diff_res",out_dir = out_path)

AEC_T1_sig = get_signature(stagepare="AEC_T1_pre_HSC", in_dir="./diff_res",out_dir = out_path)











