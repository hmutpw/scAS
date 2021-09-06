rm(list=ls())
gc()
########################################
#get teanscript infor.
########################################
gene_infor = read.table("../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T,row.names = 1,check.names=1)
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names = 1,check.names=1)
transcript_tpm = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample

stage=c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])

all_sig_iso_tab = read.table("./diff_iso_not_gene/AEC_HEC_T1/AEC_HEC_T1_pre_HSC_diff_iso_not_gene_with_all_iso.txt",sep="\t",header=T,check.names=F)
all_sig_iso_tab_not_na = na.omit(all_sig_iso_tab)
all_sig_iso_tab_stst = tapply(all_sig_iso_tab_not_na$enrich_type,all_sig_iso_tab_not_na$gene_name,function(x){
length(unique(as.character(x)))})

AEC_enrich_gene = as.character(all_sig_iso_tab_not_na[all_sig_iso_tab_not_na$enrich_type=="AEC_enrich","gene_name"])
HEC_enrich_gene = as.character(all_sig_iso_tab_not_na[all_sig_iso_tab_not_na$enrich_type=="HEC_enrich","gene_name"])
T1_enrich_gene = as.character(all_sig_iso_tab_not_na[all_sig_iso_tab_not_na$enrich_type=="T1_pre_HSC_enrich","gene_name"])

AEC_HEC_overlap = intersect(AEC_enrich_gene, HEC_enrich_gene)
AEC_HEC_overlap
AEC_T1_overlap = intersect(AEC_enrich_gene, T1_enrich_gene)
AEC_T1_overlap 
HEC_T1_overlap = intersect(HEC_enrich_gene, T1_enrich_gene)
HEC_T1_overlap

all_sig_iso_tab_stst_2 = all_sig_iso_tab_stst[which(all_sig_iso_tab_stst>=2)]
sort(all_sig_iso_tab_stst_2)

all_sig_iso_tab_stst[gene]



#################################################################################################
#plot candidates 
#################################################################################################
candi_gene = c("Sf3b1", "Ntmt1","Sec31a","E130309D02Rik")
"Rpl18a" "Sf3b1"  "Tex2" "Polr2m" "Rpl30"  "Safb2" "Tulp4"
#------get sig isoform.
gene = "E130309D02Rik"
candidates_frm = all_sig_iso_tab[all_sig_iso_tab$gene_name == gene,]
candidates_frm
candidate_gene_id = unique(as.character(candidates_frm$gene_id))
sig_iso = unique(as.character(candidates_frm[!is.na(candidates_frm$enrich_type),"transcript_id"]))[1]
not_sig_iso = setdiff(as.character(candidates_frm$transcript_id),sig_iso[1])
sig_type = c(rep("sig",length(sig_iso)),rep("not_sig",length(not_sig_iso)))
names(sig_type) = c(sig_iso, not_sig_iso)

#------get candidates iso exp
candidates_iso_exp = transcript_tpm[names(sig_type),samples]
library(data.table)
candidates_iso_exp_melt = melt(t(candidates_iso_exp))
colnames(candidates_iso_exp_melt) = c("sample","transcript_id","TPM")
candidates_iso_exp_melt$transcript_name = transcript_infor[as.character(candidates_iso_exp_melt$transcript_id),"transcript_name"]
candidates_iso_exp_melt$sig_type = factor(sig_type[as.character(candidates_iso_exp_melt$transcript_id)],
levels = rev(c("sig","not_sig")))
candidates_iso_exp_melt$stage = factor(as.character(sample2stage[as.character(candidates_iso_exp_melt$sample),"stage"]),levels = stage)

#------sort sample
sort_sig_tab = candidates_iso_exp_melt[candidates_iso_exp_melt$sig_type=="sig",]
sort_sig_iso_median = tapply(sort_sig_tab$TPM,as.character(sort_sig_tab$transcript_name),median)
sort_sig_iso_max = names(sort(sort_sig_iso_median,decreasing=T))[1]
sort_sig_tab_filter = sort_sig_tab[sort_sig_tab$transcript_name ==sort_sig_iso_max,]
sort_sig_tab_ordered = sort_sig_tab_filter[order(sort_sig_tab_filter$stage,sort_sig_tab_filter$TPM),]
sample_ordered = unique(as.character(sort_sig_tab_ordered$sample))

#------sort transcript
iso_exp_tab = unique(candidates_iso_exp_melt[,c("TPM","transcript_name","sig_type","stage")])
iso_exp_sorted = iso_exp_tab[order(iso_exp_tab$sig_type,iso_exp_tab$TPM,decreasing=T),]
iso_exp_sorted_name = unique(as.character(iso_exp_sorted$transcript_name))

#------bar plot results.
plot_frm = candidates_iso_exp_melt
plot_frm$sample = factor(plot_frm$sample, levels = sample_ordered)
plot_frm$transcript_name = factor(plot_frm$transcript_name, levels = rev(iso_exp_sorted_name))
color_fill = c("#2080DF","#f4606c")
max = 10

library(ggplot2)
text_size = 8
line_size = 0.5/1.07

p = ggplot(plot_frm, aes(x = sample, y = log2(TPM+1), fill= sig_type)) + 
geom_bar(stat = "identity",position = "stack",size=.1,width=.8, color=NA)+
scale_y_continuous(limits=c(0,max),breaks = seq(0,max,max/2))+
scale_fill_manual(values = color_fill)+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(),
axis.ticks = element_line(),
axis.title.x = element_text(size=text_size,color="black"),
axis.text.x= element_text(size=text_size,angle=90,vjust=0.5), 
axis.text.y= element_text(size=text_size,color="black",hjust=1),
legend.text = element_text(size = text_size,color="black"), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_wrap(~transcript_name,ncol=1)

p

########################################
#plot candidates gene level expression.
########################################

candi_gene_frm = gene_infor[gene_infor$gene_name %in% gene,]
candi_gene_exp = colSums(candidates_iso_exp)
candi_gene_exp_melt = data.frame(sample = names(candi_gene_exp),
gene_id = rep(candidate_gene_id,length(candi_gene_exp)),TPM = candi_gene_exp)
candi_gene_exp_melt$stage = factor(as.character(sample2stage[as.character(candi_gene_exp_melt$sample),"stage"]),levels = stage)

#------bar plot results.
plot_frm = candi_gene_exp_melt
plot_frm$sample = factor(plot_frm$sample, levels = sample_ordered)
plot_frm$gene_name = as.character(candi_gene_frm[plot_frm$gene_id,"gene_name"])
color_fill = c("#2080DF","#f4606c")
#color_fill = c("#DF2080")

max = 10

library(ggplot2)
text_size = 8
line_size = 0.5/1.07

q <- ggplot(plot_frm, aes(x = sample, y = log2(TPM+1), fill= gene_name)) + 
geom_bar(stat = "identity",position = "stack",size=.1,width=.8, color=NA)+
scale_y_continuous(limits=c(0,max),breaks = seq(0,max,max/2))+
#scale_fill_manual(values = color_fill)+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(),
axis.ticks = element_line(),
axis.title.x = element_text(size=text_size,color="black"),
axis.text.x= element_text(size=text_size,angle=90,vjust=0.5), 
axis.text.y= element_text(size=text_size,color="black",hjust=1),
legend.text = element_text(size = text_size,color="black"), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')+
facet_wrap(~gene_name,ncol=1)

q



