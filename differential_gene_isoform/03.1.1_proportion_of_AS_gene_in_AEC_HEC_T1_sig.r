rm(list=ls())
gc()
####################################################################################
#------get transcript & sample infor.
####################################################################################
library(data.table)
event_psi = read.table("../MISO_PSI/merge_PSI/merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,row.names=1,check.names=F)
event_infor = read.table("../MISO_PSI/merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1,check.names=F)

transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",row.names=1,header=T,stringsAsFactors=F)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage)= sample2stage$sample
stages = c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stages,])

########################################
#get sig isoform.
########################################
sig_path = "./DESeq2_res/transcript_level/all_res"
AEC_HEC_res = read.table(file.path(sig_path,"AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_res = na.omit(AEC_HEC_res)
AEC_HEC_diff = AEC_HEC_res[AEC_HEC_res$padj<0.05,]
AEC_HEC_sig = AEC_HEC_diff[AEC_HEC_diff$log2FoldChange>0,]

HEC_T1_sig = read.table(file.path(sig_path,"HEC_T1_sig_isoform_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_sig = read.table(file.path(sig_path,"AEC_T1_sig_isoform_res.txt"),sep="\t",header=T,row.names=1,check.names=F)

sig_iso_frm = HEC_T1_sig

sig_iso_events_in = row.names(event_infor[event_infor[,'transcript_id(in)'] %in% row.names(sig_iso_frm),])
sig_iso_events_ex = row.names(event_infor[event_infor[,'transcript_id(ex)'] %in% row.names(sig_iso_frm),])
sig_iso_events = union(sig_iso_events_in, sig_iso_events_ex)
length(sig_iso_events)

#------get PSI statistics
AEC_HEC_T1_event_psi = event_psi[,samples]

event_psi_stat = apply(AEC_HEC_T1_event_psi,2,function(x){as.character(event_infor[names(x[which(!is.na(x))]),'gene_id'])})
AS_gene_id = unique(unlist(event_psi_stat))

AEC_HEC_AS_gene = intersect(AS_gene_id,as.character(AEC_HEC_sig$gene_id))

