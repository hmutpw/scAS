rm(list=ls())
gc()
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

gene_tpm_tab = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",header=T,row.names=1,check.names=F)
gene_tpm = gene_tpm_tab[,samples]

transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]

event_psi = read.table("../MISO_PSI/merge_PSI/merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,row.names=1,check.names=F)
event_infor = read.table("../MISO_PSI/merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1,check.names=F)

#####################################
#get AS isoform overlap
#####################################
annotated_AS = intersect(row.names(event_psi), row.names(event_infor))
event_psi_annot = event_psi[annotated_AS,samples]
event_psi_in_isoform = apply(event_psi_annot,2,function(x){unlist(strsplit(as.character(event_infor[names(x[!is.na(x)]),'transcript_id(in)']),split=';'))})
event_psi_ex_isoform = apply(event_psi_annot,2,function(x){unlist(strsplit(as.character(event_infor[names(x[!is.na(x)]),'transcript_id(ex)']),split=';'))})

#------get expressed isoform (TPM>1).
getnum <- function(x ,cutoff = 1, iso_info = transcript_infor){
iso_id = names(x[which(x > cutoff)])
iso_exp_tab = iso_info[iso_id,]
gene_stat = table(iso_exp_tab[iso_id,"gene_id"])
iso_exp_tab$iso_num = gene_stat[as.character(iso_exp_tab$gene_id)]
two_iso_unique = row.names(iso_exp_tab[iso_exp_tab$iso_num>1,])
return(two_iso_unique)
}
multi_iso_name = apply(transcript_tpm,2,FUN=getnum)
multi_iso_num = sapply(multi_iso_name,length)

AS_multi_iso_num=c()
for(i in names(multi_iso_num)){
AS_in_iso = event_psi_in_isoform[[i]]
AS_ex_iso = event_psi_ex_isoform[[i]]
multi_iso_id = multi_iso_name[[i]]

AS_in_multi_iso = intersect(AS_in_iso, multi_iso_id)
AS_ex_multi_iso = intersect(AS_ex_iso, multi_iso_id)
AS_in_ex_multi_iso = union(AS_in_multi_iso, AS_ex_multi_iso)
AS_in_ex_multi_iso_num = length(AS_in_ex_multi_iso)
names(AS_in_ex_multi_iso_num) = i
AS_multi_iso_num=c(AS_multi_iso_num, AS_in_ex_multi_iso_num)
}


multi_iso_tab = data.frame(multi_iso = multi_iso_num,AS_iso = AS_multi_iso_num[names(multi_iso_num)])
multi_iso_tab$not_AS_iso = apply(multi_iso_tab,1,function(x){x[1]-x[2]})

multi_iso_tab_for_sort = multi_iso_tab
multi_iso_tab_for_sort$stage = factor(sample2stage[row.names(multi_iso_tab_for_sort),"stage"],levels = rev(stage))
multi_iso_tab_for_sort_order = multi_iso_tab_for_sort[order(multi_iso_tab_for_sort$stage,multi_iso_tab_for_sort$AS_iso,decreasing=T),]
ordered_sample = row.names(multi_iso_tab_for_sort_order)

library(data.table)
multi_iso_melt = melt(t(multi_iso_tab[,c("AS_iso","not_AS_iso")]))
colnames(multi_iso_melt) = c("type","sample","num")
multi_iso_melt$sample = factor(multi_iso_melt$sample,levels = ordered_sample)
multi_iso_melt$type = factor(multi_iso_melt$type,levels = rev(c("AS_iso","not_AS_iso")))

plot_frm = multi_iso_melt
library(ggplot2)

text_size = 8
line_size = 0.5
dodge = position_dodge(0)

p = ggplot(plot_frm,aes(x = sample, y=num, fill=type))+
geom_bar(stat="identity",color = NA,alpha=0.75)+
scale_y_continuous(limits=c(0,11000),breaks=seq(0,10000,by=2500))+
scale_fill_manual(values=c("#2020DF","#DF20DF"))+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=0.5*line_size, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45,hjust=1,vjust=1,color="black"), 
axis.text.y= element_text(size=text_size, color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p



