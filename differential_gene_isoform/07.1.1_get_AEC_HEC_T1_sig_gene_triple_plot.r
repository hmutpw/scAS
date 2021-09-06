rm(list=ls())
gc()
########################################
#get teanscript infor.
########################################
gene_infor = read.table("../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",row.names=1,header=T,stringsAsFactors=F)
gene_tpm = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample

stage=c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])

gene_tpm_median = t(apply(gene_tpm[,samples],1,function(x){tapply(x,as.character(sample2stage[names(x),"stage"]),median)}))

gene_dir = "./diff_res/gene_level/pairwise_test/all_res"
AEC_HEC_res = read.table(file.path(gene_dir,"AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1)
AEC_HEC_high = AEC_HEC_res[apply(AEC_HEC_res[,grep("_median_TPM$",colnames(AEC_HEC_res))],1,function(x){length(x[which(x>1)])>0}),]
AEC_HEC_diff = AEC_HEC_high[AEC_HEC_high$FDR<0.05,]
AEC_HEC_diff_up = AEC_HEC_diff[AEC_HEC_diff$logFC>0,]
AEC_HEC_diff_down = AEC_HEC_diff[AEC_HEC_diff$logFC<0,]

HEC_T1_res = read.table(file.path(gene_dir,"HEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1)
HEC_T1_high = HEC_T1_res[apply(HEC_T1_res[,grep("_median_TPM$",colnames(HEC_T1_res))],1,function(x){length(x[which(x>1)])>0}),]
HEC_T1_diff = HEC_T1_high[HEC_T1_high$FDR<0.05,]
HEC_T1_diff_up = HEC_T1_diff[HEC_T1_diff$logFC>0,]

#------
HEC_sig = row.names(AEC_HEC_diff_up)
HEC_enrich_gene = rep("HEC signature gene",length(HEC_sig))
names(HEC_enrich_gene) = HEC_sig

T1_sig = setdiff(row.names(HEC_T1_diff_up), HEC_sig)
T1_enrich_gene = rep("T1 pre-HSC signature gene",length(T1_sig))
names(T1_enrich_gene) = T1_sig

AEC_sig = setdiff(row.names(AEC_HEC_diff_down),union(HEC_sig, T1_sig))
AEC_enrich_gene = rep("AEC upregulated",length(AEC_sig))
names(AEC_enrich_gene) = AEC_sig

all_sig_gene = c(AEC_enrich_gene, HEC_enrich_gene, T1_enrich_gene)

all_sig_gene_frm = data.frame(gene_id = names(all_sig_gene),gene_name = AEC_HEC_res[names(all_sig_gene),"gene_name"],
gene_tpm_median[names(all_sig_gene),],enrich_type = all_sig_gene)

write.table(all_sig_gene_frm,"./diff_res/gene_level/pairwise_test/diff_res/AEC_HEC_T1_sig_gene_type_table.txt",sep="\t",quote=F,row.names=F)


candi_gene = c("Runx1","Dnmt3b")
candi_plot_frm = all_sig_gene_frm[all_sig_gene_frm$gene_name %in% candi_gene,]
candi_plot_frm$enrich_type = candi_plot_frm$gene_name
candi_plot_frm_all = all_sig_gene_frm[setdiff(row.names(all_sig_gene_frm),row.names(candi_plot_frm)),]
candi_plot_frm_all$gene_name = rep("",nrow(candi_plot_frm_all))
candi_plot_frm_all_new = rbind(candi_plot_frm_all,candi_plot_frm)
plot_tab = candi_plot_frm_all_new

library("ggtern")
text_size = 8
line_size = 0.5/1.07

p1<-ggtern(data=plot_tab,aes(x=HEC,y=AEC,z=T1_pre_HSC, color= enrich_type,label = gene_name))+
geom_point(alpha=1,size=2,shape=16)+
geom_text(color="black",size=2)+
scale_color_manual(values = c("#1565c0","#9e9d24","#d84315",rep("purple",10)))+
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
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p1

