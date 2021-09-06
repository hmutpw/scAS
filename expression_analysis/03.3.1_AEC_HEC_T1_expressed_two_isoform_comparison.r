rm(list=ls())
gc()
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
gene_infor = read.table("../ref_genome/gencode.vM22.gene.infor.rm.version.tsv",sep="\t",header=T)
row.names(gene_infor) = gene_infor$gene_id
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

transcript_tpm_tab = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",header=T,row.names=1,check.names=F)
transcript_tpm = transcript_tpm_tab[,samples]
gene_num_mat = apply(transcript_tpm,2,function(x){table(transcript_infor[names(which(x>1)),'gene_id'])})

###########################################
#------get multi isoform gene/isoform num
###########################################
gene_num_type = gene_num_mat
#if the expressed isoform number >1,asign it to 2.
gene_num_type = ifelse(gene_num_type>1,2,gene_num_type)
gene_num_type_filter = gene_num_type[rowSums(gene_num_type)>0,]

# if a gene in this stage expressed more than 2 isoforms (TPM>1) in at least 5
# cells, we consider it as a gene with multiple type in this stage. 
assign_type <- function(x){
gene_num_type = NULL
stst_mat = rep(0,3)
names(stst_mat) = c(0:2)
stst_mat[names(x)] = x
if(stst_mat['2']>=5){
gene_num_type = "multiple"
}else if(sum(stst_mat['1'],stst_mat['2'])>=5){
if(stst_mat['1']>stst_mat['2']){
gene_num_type = "single"
}else{
gene_num_type = "multiple"
}
}else{
gene_num_type = 'no expression'
}
return(gene_num_type)
}

#------merge gene type for each stage.
multi_iso_gene_stage = data.frame(gene_id = row.names(gene_num_type_filter))
for(i in stage){
stage_samp = as.character(sample2stage[sample2stage$stage %in% i,"sample"])
stage_iso_gene_mat = gene_num_type_filter[,stage_samp]
stage_iso_gene_stat = apply(stage_iso_gene_mat,1,table)
stage_iso_gene_type = sapply(stage_iso_gene_stat,assign_type)
multi_iso_gene_frm = data.frame(gene_id = names(stage_iso_gene_type),type = stage_iso_gene_type)
colnames(multi_iso_gene_frm)[2] = i
multi_iso_gene_stage = merge(multi_iso_gene_stage,multi_iso_gene_frm,by="gene_id",all=T,sorted=F)
}
#write.table(multi_iso_gene_stage,"./AEC_HEC_T1_multi_iso/7_stage_all_3_gene_type_table.txt",sep="\t",row.names=F,quote=F)

multi_iso_gene_stat = multi_iso_gene_stage[,-1]
row.names(multi_iso_gene_stat) = multi_iso_gene_stage[gene_infor[],1]

###########################
#get two types of multi isoform table.
###########################
AEC_not_multi_gene_id = row.names(multi_iso_gene_stat[multi_iso_gene_stat$AEC!="multiple",])
HEC_not_multi_gene_id = row.names(multi_iso_gene_stat[multi_iso_gene_stat$HEC!="multiple",])
HEC_multi_gene_id = row.names(multi_iso_gene_stat[multi_iso_gene_stat$HEC=="multiple",])
T1_multi_gene_id = row.names(multi_iso_gene_stat[multi_iso_gene_stat$T1_pre_HSC=="multiple",])
AEC_not_HEC_T1_multi_gene = Reduce(intersect,list(AEC_not_multi_gene_id, HEC_multi_gene_id,T1_multi_gene_id))
AEC_HEC_not_T1_multi_gene = Reduce(intersect,list(AEC_not_multi_gene_id, HEC_not_multi_gene_id,T1_multi_gene_id))
T1_multi_gene_type = c(rep("HEC_T1_multi",length(AEC_not_HEC_T1_multi_gene)),
rep("T1_only_multi",length(AEC_HEC_not_T1_multi_gene)))
names(T1_multi_gene_type) = c(AEC_not_HEC_T1_multi_gene, AEC_HEC_not_T1_multi_gene)

HEC_T1_multi_gene_id = union(AEC_not_HEC_T1_multi_gene, AEC_HEC_not_T1_multi_gene)
HEC_T1_multi_tab = multi_iso_gene_stat[HEC_T1_multi_gene_id,]
HEC_T1_multi_tab_out = data.frame(gene_infor[row.names(HEC_T1_multi_tab),c("gene_id","gene_name")],HEC_T1_multi_tab)
write.table(HEC_T1_multi_tab_out,"./AEC_HEC_T1_multi_iso/7_stage_single_multi_gene.txt",sep="\t",row.names=F,quote=F)


#HEC_T1_multi_tab$T1_pre_HSC = T1_multi_gene_type[row.names(HEC_T1_multi_tab)]

multi_iso_gene_stat_table = as.data.frame(table(as.data.frame(HEC_T1_multi_tab)))
sort_name = c("HEC_T1_multi", "T1_only_multi", "multiple","single","no expression")

data_edges = multi_iso_gene_stat_table
edges_filter = data_edges[data_edges$Freq>0,]
library(data.table)
data = melt(edges_filter)

library(ggalluvial)
library(ggplot2)
library(ggforce)
data = edges_filter
data <- gather_set_data(data, 1:length(stage))
data$x = factor(data$x,levels=stage)
data$y = factor(data$y,levels=sort_name)
data$HEC = factor(data$HEC,levels=sort_name)

stage_color = c("#e53935","#1e88e5","#757575")
p = ggplot(data, aes(x = x, id = id, split = y, value = Freq)) +
geom_parallel_sets(aes(fill = Adult_HSC), alpha = 0.5, axis.width = 0.1) +
geom_parallel_sets_axes(axis.width = 0.1) +
geom_parallel_sets_labels(colour = 'white')
+
scale_color_manual(values=stage_color)+
scale_fill_manual(values=stage_color)
p

p = ggplot(data, aes(axis1 = AEC, axis2 = HEC, axis3 = T1_pre_HSC,
axis4 =  T2_pre_HSC, axis5 = E12, axis6 = E14, axis7 = Adult_HSC,weight = Freq,fill=AEC)) +
  geom_flow() +
  geom_stratum(alpha = .5) 
p

data_new = multi_iso_gene_stat
T1_multi = data_new[data_new$AEC!="multiple",]
HEC_T1_HSC_multi = T1_multi[apply(T1_multi[,-1],1,function(x){length(x[which(x=="multiple")])==length(x)}),]
HEC_T1_HSC_multi_type = rep("HEC_T1_HSC",nrow(HEC_T1_HSC_multi))
names(HEC_T1_HSC_multi_type) = row.names(HEC_T1_HSC_multi)
T1_only = T1_multi[T1_multi$HEC!="multiple",]
T1_only_multi = T1_only[apply(T1_only[,-c(1:2)],1,function(x){length(x[which(x=="multiple")])==length(x)}),]
T1_only_multi_type = rep("T1_HSC",nrow(T1_only_multi))
names(T1_only_multi_type) = row.names(T1_only_multi)

all_type = c(HEC_T1_HSC_multi_type, T1_only_multi_type)

HEC_T1_multi_tab_new = HEC_T1_multi_tab
HEC_T1_multi_tab_new_1 = HEC_T1_multi_tab_new[names(all_type),]
HEC_T1_multi_tab_new_1$Adult_HSC = all_type[row.names(HEC_T1_multi_tab_new_1)]
HEC_T1_multi_tab_new_2 = HEC_T1_multi_tab_new[setdiff(row.names(HEC_T1_multi_tab_new),names(all_type)),]

HEC_T1_multi_tab_new_tab = rbind(HEC_T1_multi_tab_new_2, HEC_T1_multi_tab_new_1)

multi_iso_gene_stat_table = as.data.frame(table(as.data.frame(HEC_T1_multi_tab_new_tab)))
sort_name = c("HEC_T1_HSC", "T1_HSC", "multiple","single","not_expressed")

data_edges = multi_iso_gene_stat_table
edges_filter = data_edges[data_edges$Freq>0,]


library(ggalluvial)
library(ggplot2)
library(ggforce)
data = edges_filter
data <- gather_set_data(data, 1:length(stage))
data$x = factor(data$x,levels=stage)
data$y = factor(data$y,levels=sort_name)
data$Adult_HSC = factor(data$HEC,levels=sort_name)

stage_color = c("#e53935","#1e88e5","#757575")
p = ggplot(data, aes(x = x, id = id, split = y, value = Freq)) +
geom_parallel_sets(aes(fill = T1_pre_HSC), alpha = 0.5, axis.width = 0.1) +
geom_parallel_sets_axes(axis.width = 0.1) +
geom_parallel_sets_labels(colour = 'white')+
scale_color_manual(values=stage_color)+
scale_fill_manual(values=stage_color)
p



