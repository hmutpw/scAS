rm(list=ls())
gc()
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])


transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]
transcript_exp_num = apply(transcript_tpm,2,function(x){length(x[which(x>1)])})
transcript_num_frm = data.frame(sample = names(transcript_exp_num),
number = transcript_exp_num,type=rep("isoform",length(transcript_exp_num)))

gene_tpm_tab = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",header=T,row.names=1,check.names=F)
gene_tpm = gene_tpm_tab[,samples]
gene_exp_num = apply(gene_tpm,2,function(x){length(x[which(x>1)])})
write.table(gene_exp_num,"expressed_gene_number.txt",sep="\t",quote=F)
gene_num_frm = data.frame(sample = names(gene_exp_num),
number = gene_exp_num,type=rep("gene",length(gene_exp_num)))

plot_frm = rbind(transcript_num_frm, gene_num_frm)
plot_frm$stage = sample2stage[as.character(plot_frm$sample),"stage"]
plot_frm$stage = factor(plot_frm$stage,levels=stage)
plot_frm$type = factor(plot_frm$type,levels=c("gene","isoform"))
plot_frm$x = apply(plot_frm,1,function(x){paste(x["type"],x["stage"],sep="_")})
levels_x = paste(c("gene","isoform"),rep(stage,each=2),sep="_")
plot_frm$x = factor(plot_frm$x,levels = levels_x)

iso_num_median = tapply(plot_frm$number,plot_frm$type,median)

text_size = 8
line_size = 0.5/1.07

library(ggplot2)
p <- ggplot(plot_frm,mapping = aes(x=x,y=number/1000,fill=type))+
geom_hline(yintercept = c(iso_num_median/1000))+
geom_violin(position = "dodge", na.rm = TRUE,scale = "width",color=NA,width = 0.9)+
geom_jitter(width=.2,size= .3,shape=16,color="black")+
scale_fill_manual(values=c("#3288BD","#F781BF"))+
#scale_color_manual(values=stage2color[sixStage,"figurecolor"])+
scale_y_continuous(limits = c(0,17), breaks = seq(0,16,by=4))+
labs(x="",y='number of expressed gene/isoform in each single cell (X10^3)')+
theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p


isoform_tab = plot_frm[plot_frm$type=="isoform",]




AEC_frm = isoform_tab[isoform_tab$stage=='AEC',"number"]
HEC_frm = isoform_tab[isoform_tab$stage=='HEC',"number"]
T1_frm = isoform_tab[isoform_tab$stage=='T1_pre_HSC',"number"]
T2_frm = isoform_tab[isoform_tab$stage=='T2_pre_HSC',"number"]
AEC_HEC = wilcox.test(HEC_frm,AEC_frm, alternative ="greater")
HEC_T1 = wilcox.test(T1_frm, HEC_frm, alternative ="greater")
T1_T2 = wilcox.test(T1_frm, T2_frm, alternative ="greater")


