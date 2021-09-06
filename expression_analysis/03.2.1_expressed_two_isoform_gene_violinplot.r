rm(list=ls())
gc()
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]
gene_num_mat = apply(transcript_tpm,2,function(x){table(transcript_infor[names(which(x>1)),'gene_id'])})
gene_num = apply(gene_num_mat,2,function(x){c(table(x)[1:2],'>1'=sum(table(x)[-c(1:2)]))})

###########################################
#------get multi isoform gene/isoform num
###########################################
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
multi_iso_num_frm = data.frame(sample = names(multi_iso_num), num = multi_iso_num)
multi_iso_num_frm$stage = sample2stage[as.character(multi_iso_num_frm$sample),"stage"]
multi_iso_num_frm$type = rep("isoform",nrow(multi_iso_num_frm))
#------gene level.
multi_iso_gene_num = gene_num['>1',]
multi_iso_gene_frm = data.frame(sample = names(multi_iso_gene_num), num = multi_iso_gene_num)
multi_iso_gene_frm$stage = sample2stage[as.character(multi_iso_gene_frm$sample),"stage"]
multi_iso_gene_frm$type = rep("gene",nrow(multi_iso_gene_frm))

multi_iso_num_plot = rbind(multi_iso_gene_frm,multi_iso_num_frm)
multi_iso_num_plot$stage = factor(multi_iso_num_plot$stage,levels = stage)
multi_iso_num_plot$type = factor(multi_iso_num_plot$type,levels = c("gene","isoform"))
multi_iso_num_plot = multi_iso_num_plot[order(multi_iso_num_plot$stage),]
multi_iso_num_plot$id = apply(multi_iso_num_plot,1,function(x){paste(x['type'],x['stage'],sep="_")})
multi_iso_num_plot$id = factor(multi_iso_num_plot$id,levels = (unique(as.character(multi_iso_num_plot$id))))

multi_iso_gene_num_frm = multi_iso_num_frm
multi_iso_gene_num_frm$gene_num = multi_iso_gene_frm[row.names(multi_iso_gene_num_frm),"num"]
multi_iso_gene_num_frm$ratio = apply(multi_iso_gene_num_frm,1,function(x){as.numeric(x["num"])/as.numeric(x["gene_num"])})

multi_iso_num_median = tapply(multi_iso_num_plot$num,multi_iso_num_plot$type,median)

plot_frm = multi_iso_gene_num_frm
plot_frm$stage = factor(plot_frm$stage,levels = stage)
sort_sample = as.character(plot_frm[order(plot_frm$stage),'sample'])
plot_frm$sample = factor(plot_frm$sample,levels = sort_sample)


text_size = 8
line_size = 0.5/1.07

library(ggplot2)
p <- ggplot(plot_frm,mapping = aes(x=sample,y=ratio,fill = stage))+
#geom_hline(yintercept = c(multi_iso_num_median/1000))+
#geom_violin(position = "dodge", na.rm = TRUE,scale = "width",color=NA,width = 0.9)+
#geom_jitter(width=0.2,size= .3,shape=16,color="black")+
geom_bar(stat = "identity")+
#scale_fill_manual(values=c("#3288BD","#F781BF"))+
#scale_color_manual(values=stage2color[sixStage,"figurecolor"])+
#scale_y_continuous(limits = c(0,12), breaks = seq(0,12,by=3))+
labs(x="",y='isoform/gene')+
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






















