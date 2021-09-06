rm(list=ls())
gc()

transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]
gene_tpm_tab = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",header=T,row.names=1,check.names=F)
gene_tpm = gene_tpm_tab[,samples]

gene_num_mat = apply(transcript_tpm,2,function(x){table(transcript_infor[names(which(x>1)),'gene_id'])})
gene_num_1 = sapply(gene_num_mat,function(x){length(x[which(x==1)])})
gene_num_2 = sapply(gene_num_mat,function(x){length(x[which(x>1)])})
gene_num = data.frame(sample = names(gene_num_1), gene_num_1,gene_num_2)
gene_num$ratio = apply(gene_num,1,function(x){as.numeric(x[3])/(as.numeric(x[2])+as.numeric(x[3]))})
gene_num$stage = factor(sample2stage[as.character(gene_num$sample),"stage"],levels = stage)
gene_num$sample = factor(gene_num$sample,levels = samples)

#-----proportion of expressed isoform with >=2 isoform gene ratio.
library(data.table)
plot_tab = gene_num

median_ratio = median(plot_tab$ratio)

stage2color = read.table("../sample_infor/stage2color.txt",sep="\t",header=T,comment.char = "",stringsAsFactors = F,row.names=1)

text_size = 8
line_size = 0.5/1.07

library(ggplot2)
p<-ggplot(plot_tab, aes(x = sample ,y = ratio, color = stage, fill = stage))+
geom_point(size=1, shape=16)+
geom_hline(yintercept = c(median_ratio,0.5))+
scale_color_manual(values=stage2color[stage,"figurecolor"])+
scale_fill_manual(values=stage2color[stage,"figurecolor"])+
scale_y_continuous(limits = c(0,0.5), breaks = seq(0,0.5,by=.1))+
labs(x="",y='isoform proportaion')+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),
legend.position='right')

p


T1_frm = plot_tab[plot_tab$stage=="T1_pre_HSC",]
T1_frm_1 = T1_frm[T1_frm$type=="1","ratio"]
T1_frm_more_1 = T1_frm[T1_frm$type==">1","ratio"]

other_stages = stage[which(stage!="T1_pre_HSC")]
other_stages_frm = plot_tab[plot_tab$stage %in% other_stages,]
other_stages_frm_1 = other_stages_frm[other_stages_frm$type=="1","ratio"]
other_stages_frm_more_1 = other_stages_frm[other_stages_frm$type==">1","ratio"]

T1_others_1_test = wilcox.test(T1_frm_1,other_stages_frm_1,exact=T,alternative="less")$p.value
T1_others_more_1_test = wilcox.test(T1_frm_more_1,other_stages_frm_more_1,exact=T,alternative="greater")$p.value

P_values = c()
Aver_ratio <- c()
for(i in other_stages){
stage_frm = plot_tab[plot_tab$stage==i,]
#stage_frm_more_1 = stage_frm[stage_frm$type==">1","ratio"]
average_ratio <- median(T1_frm$ratio)-median(stage_frm$ratio)
names(average_ratio) <- i
p = wilcox.test(T1_frm$ratio, stage_frm$ratio, exact=T,alternative="greater")$p.value
names(p)=i
P_values = c(P_values, p)
Aver_ratio <- c(Aver_ratio, average_ratio)
}
P_values
Aver_ratio


write.table(P_values,"./figure/isoform_ratio_T1_other_stages_wc_test_p_value.txt",sep="\t",col.names=F,quote=F)


median(plot_tab$gene_num_2)

median(T1_frm$gene_num_2)-median(other_stages_frm$gene_num_2)











