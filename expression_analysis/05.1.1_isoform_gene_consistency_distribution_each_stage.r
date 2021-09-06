rm(list=ls())
gc()
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]

stage_gene_multi = read.table("./AEC_HEC_T1_multi_iso/7_stage_all_3_gene_type_table.txt",sep="\t",header=T,row.names=1)
stage_gene_multi_gene = stage_gene_multi[apply(stage_gene_multi,1,function(x){"multiple" %in% x}),]

######################################################################################################
#expression level
######################################################################################################

get_cor <- function(stage, iso_exp = transcript_tpm, isoformInfor = transcript_infor, samp2Stage=sample2stage){

#------get samples of stage.
samples = as.character(samp2Stage[samp2Stage$stage==stage,"sample"])
samp_num = length(samples)

#------get isoform/gene expression.
isoform_exp = intersect(row.names(iso_exp),row.names(isoformInfor))
stage_iso_exp = iso_exp[isoform_exp,samples]
#------isoform TPM>1 in at least 4 cells each stage were remined for cor analysis.
getnum <- function(x ,cutoff = 1, iso_info = isoformInfor){
iso_id = names(x[which(x > cutoff)])
iso_exp_tab = iso_info[iso_id,]
gene_stat = table(iso_exp_tab[iso_id,"gene_id"])
iso_exp_tab$iso_num = gene_stat[as.character(iso_exp_tab$gene_id)]
two_iso_unique = row.names(iso_exp_tab[iso_exp_tab$iso_num>1,])
return(two_iso_unique)
}
multi_iso_name = apply(stage_iso_exp,2,FUN=getnum)
multi_iso_id = table(unlist(multi_iso_name))
stage_multi_iso = names(multi_iso_id[which(multi_iso_id>=5)])
stage_multi_iso_exp = stage_iso_exp[stage_multi_iso,]
#------get gene exp.
stage_gene_exp = apply(stage_multi_iso_exp,2,function(x)(tapply(x,as.character(isoformInfor[names(x),"gene_id"]),sum)))
#------get gene isoform consistency.
merge_iso_gene_frm = data.frame(stage_multi_iso_exp,
stage_gene_exp[as.character(isoformInfor[row.names(stage_multi_iso_exp),"gene_id"]),],check.names = F)
stage_iso_gene_cor = apply(merge_iso_gene_frm,1,
function(x){cor(as.numeric(x[1:samp_num]),as.numeric(x[(samp_num+1):ncol(merge_iso_gene_frm)]))})
stage_gene_median_cor = tapply(stage_iso_gene_cor,as.character(isoformInfor[names(stage_iso_gene_cor),"gene_id"]),median)

out_frm = data.frame(gene_id = names(stage_gene_median_cor),cor = stage_gene_median_cor,
stage = rep(stage,length(stage_gene_median_cor)),check.names = F)
return(out_frm)
}

cor_frm = c()
for(i in stage){
cor_tab = get_cor(stage=i)
cor_frm = rbind(cor_frm, cor_tab)
}

############################################################
stage2color = read.table("../sample_infor/stage2color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)

allStageCol = stage2color[stage,"figurecolor"]
library(ggplot2)
library(splines)
library(plyr)

plot_frm = ddply(cor_frm, .(stage), transform, ecdf = ecdf(cor)(cor) )

text_size = 8
line_size = 0.5
p = ggplot(plot_frm,aes(x = cor, y=ecdf, color = stage))+
stat_ecdf(size = line_size)+
labs(x = 'Pearson Correlation Coefficient', y = 'Cumulative Distribution')+
#xlim(0.5,1.0)+ylim(0.5,0.75)+
#xlim(0,1.0)+
scale_color_manual(values=allStageCol)+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=0.5*line_size, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, hjust=0.5,vjust=1,color="black"), 
axis.text.y= element_text(size=text_size, color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p


#ks.test

T1_pre_HSC_cor = cor_frm[cor_frm$stage=="T1_pre_HSC","cor"]
T2_pre_HSC_cor = cor_frm[cor_frm$stage=="T2_pre_HSC","cor"]
AEC_cor = cor_frm[cor_frm$stage=="AEC","cor"]
HEC_cor = cor_frm[cor_frm$stage=="HEC","cor"]


ks.test(T1_pre_HSC_cor, T2_pre_HSC_cor)$p.value




