rm(list=ls())
gc()

gene_infor <- read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
gene_tpm <- read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",sep="\t",header=T,row.names=1)
expressed_gene_infor <- gene_infor[row.names(gene_tpm),]

single_iso_gene <- names(table(expressed_gene_infor$gene_id)[table(expressed_gene_infor$gene_id)==1])
multi_iso_gene <- setdiff(expressed_gene_infor$gene_id,single_iso_gene)

iso_path = "./diff_res/transcript_level/pairwise_test/all_res"
gene_path = "./diff_res/gene_level/pairwise_test/all_res"

stagepare = "AEC_HEC"
stagepare = "HEC_T1_pre_HSC"
exp_cutoff = 1
#------isoform level.
iso_files = file.path(iso_path, paste(stagepare,"_test_all_res.txt",sep=""))
iso_tab = read.table(iso_files,sep="\t",header=T,row.names=1,stringsAsFactors=F)
iso_high = iso_tab[apply(iso_tab[,grep("_median_TPM$",colnames(iso_tab))],1,function(x){length(x[which(x>exp_cutoff)])>0}),]
diff_iso = iso_high[iso_high$FDR<0.05,]
diff_iso = na.omit(diff_iso)
diff_iso_up = diff_iso[diff_iso$logFC>0,]
diff_iso_down = diff_iso[diff_iso$logFC<0,]

#------gene level.
gene_files = file.path(gene_path, paste(stagepare,"_test_all_res.txt",sep=""))
gene_tab = read.table(gene_files,sep="\t",header=T,row.names=1,stringsAsFactors=F)
gene_high = gene_tab[apply(gene_tab[,grep("_median_TPM$",colnames(gene_tab))],1,function(x){length(x[which(x>exp_cutoff)])>0}),]
diff_gene = gene_high[gene_high$FDR<0.05,]
diff_gene = na.omit(diff_gene)
diff_gene_up = diff_gene[diff_gene$logFC>0,]
diff_gene_down = diff_gene[diff_gene$logFC<0,]
not_diff_gene = gene_tab[gene_tab$FDR>0.1,]

#------get sig iso or sig gene
pos_up_iso_gene_id <- unique(as.character(diff_iso_up$gene_id))
sig_gene_id <- row.names(diff_gene_up)

sig_iso_not_gene_id <- intersect(pos_up_iso_gene_id,union(row.names(not_diff_gene),row.names(diff_gene_down)))
length(intersect(sig_iso_not_gene_id, multi_iso_gene))
sig_iso_not_gene_iso_id <- row.names(diff_iso_up[diff_iso_up$gene_id %in% sig_iso_not_gene_id,])
sig = length(sig_iso_not_gene_iso_id)
others <- length(setdiff(row.names(diff_iso_up), sig_iso_not_gene_iso_id))
sig/sum(sig, others)
pie(c(length(sig_iso_not_gene_iso_id), length(setdiff(row.names(diff_iso_up), sig_iso_not_gene_iso_id))))

sig_iso_gene_id <- intersect(pos_up_iso_gene_id,row.names(diff_gene_up))
sig_iso_gene_iso_id <- row.names(diff_iso_up[diff_iso_up$gene_id %in% sig_iso_gene_id,])
length(sig_iso_gene_iso_id)
sig/sum(length(sig_iso_gene_iso_id),sig)

sig_iso_gene_id_multi <- intersect(sig_iso_gene_id, multi_iso_gene)
sig_iso_gene_iso_id_multi <- row.names(diff_iso_up[diff_iso_up$gene_id %in% sig_iso_gene_id_multi,])
length(sig_iso_not_gene_iso_id)/sum(length(sig_iso_not_gene_iso_id),length(sig_iso_gene_iso_id_multi))

sig_iso_gene_id_single <- intersect(sig_iso_gene_id, single_iso_gene)
sig_iso_gene_iso_id_single <- row.names(diff_iso_up[diff_iso_up$gene_id %in% sig_iso_gene_id_single,])
length(sig_iso_not_gene_iso_id)/sum(length(sig_iso_not_gene_iso_id),length(sig_iso_gene_iso_id_multi))



pie(c(555,2833-555))
c(555)/sum(2833)




