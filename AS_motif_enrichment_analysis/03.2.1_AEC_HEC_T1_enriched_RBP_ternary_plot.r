rm(list=ls())
gc()
#install.packages("ggtern")
library("ggtern")
gene_infor = read.table("../../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
gene_tpm = read.table("../../expression_data/merged_exp/7_stage_gene_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample

stage=c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])


#------get differential test result.
test_path = "../../DEAnalysis/DESeq2_res/gene_level/all_res"
AEC_HEC_res = read.table(file.path(test_path,"AEC_HEC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_res = read.table(file.path(test_path,"HEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_res = read.table(file.path(test_path,"AEC_T1_pre_HSC_test_all_res.txt"),sep="\t",header=T,row.names=1,check.names=F)

diff_path = "../../DEAnalysis/DESeq2_res/gene_level/diff_res"
AEC_HEC_diff = read.table(file.path(diff_path,"AEC_HEC_test_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
HEC_T1_diff = read.table(file.path(diff_path,"HEC_T1_pre_HSC_test_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_T1_diff = read.table(file.path(diff_path,"AEC_T1_pre_HSC_test_diff_res.txt"),sep="\t",header=T,row.names=1,check.names=F)
AEC_HEC_T1_diff_gene = Reduce(union,list(row.names(AEC_HEC_diff), row.names(HEC_T1_diff),row.names(AEC_T1_diff)))

#------get RBP genes.
all_rMAPS_id = read.table("./rMAPS_enrichment/96_motif_RBP/rMAPS_motif_2_RBP_ID.txt",sep="\t",header=T,check.names=F)
row.names(all_rMAPS_id) = all_rMAPS_id$motif_id
RBP_list = read.table("./rMAPS_enrichment/EHT_enriched_RBPs.txt",sep="\t",header=T)
row.names(RBP_list) = RBP_list$gene_id

diff_RBP_tpm = gene_tpm[row.names(RBP_list),samples]
diff_RBP_total_median_tpm = apply(diff_RBP_tpm,1,median)
diff_RBP_tpm_median = t(apply(diff_RBP_tpm,1,function(x){tapply(x,as.character(sample2stage[names(x),"stage"]),median)}))
diff_RBP_tpm_median = diff_RBP_tpm_median[rowSums(diff_RBP_tpm_median)>0,]

plot_data = data.frame(RBP_list[row.names(diff_RBP_tpm_median),],
diff_RBP_tpm_median,
median_exp = diff_RBP_total_median_tpm[row.names(diff_RBP_tpm_median)])

plot_data$locs = factor(plot_data$locs,levels = c("exon","intron"))


p1<-ggtern(data=plot_data,aes(x=HEC,y=AEC,z=T1_pre_HSC,shape = locs))+
geom_point(alpha=.8,size=2)+
#geom_text(aes(label=gene_name),color="black",size=8)+
scale_color_manual(values = c("#DF2080","#20DF20","#2020DF"))+
scale_shape_manual(values = c(16,17,15))

p1

write.table(plot_data,"AEC_HEC_T1_diff_RBPs.txt",sep="\t",quote=F,row.names=F)


df = data.frame(x = runif(50),
 y = runif(50),
 z = runif(50),
 f = seq(1:50),
 Value = runif(50,1,10),
 Group = as.factor(round(runif(50,1,2))))
labs <- labs(x = "X", y = "Y", z = "Z", title = "Title")

P_Tern <- ggtern(data = df,aes(x, y, z)) + 
 geom_point(aes(color = Group)) + 
 labs()

