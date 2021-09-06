rm(list=ls())
gc()
########################################
#get teanscript infor.
########################################
transcript_tpm = read.table("../expression_data/merged_exp/7_stage_transcript_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
gene_tpm = read.table("../expression_data/merged_exp/7_stage_gene_TPM.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample

stage=c("AEC","HEC","T1_pre_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stage,])

in_dir = "./diff_iso_and_gene/AEC_HEC_T1"
AEC_HEC_diff_iso_and_gene_tab = read.table(file.path(in_dir,"AEC_HEC_diff_iso_and_gene.txt"),sep="\t",header=T,check.names=F)
row.names(AEC_HEC_diff_iso_and_gene_tab) = AEC_HEC_diff_iso_and_gene_tab$transcript_id
AEC_HEC_diff_iso_and_gene_tab = na.omit(AEC_HEC_diff_iso_and_gene_tab)
AEC_enriched = AEC_HEC_diff_iso_and_gene_tab[AEC_HEC_diff_iso_and_gene_tab$enrich_type=="AEC_enrich",
c("transcript_id","transcript_name","gene_id","gene_name","enrich_type")]
HEC_enriched = AEC_HEC_diff_iso_and_gene_tab[AEC_HEC_diff_iso_and_gene_tab$enrich_type=="HEC_enrich",
c("transcript_id","transcript_name","gene_id","gene_name","enrich_type")]
HEC_T1_diff_iso_and_gene_tab = read.table(file.path(in_dir,"HEC_T1_pre_HSC_diff_iso_and_gene.txt"),sep="\t",header=T,check.names=F)
row.names(HEC_T1_diff_iso_and_gene_tab) = HEC_T1_diff_iso_and_gene_tab$transcript_id
HEC_T1_diff_iso_and_gene_tab = na.omit(HEC_T1_diff_iso_and_gene_tab)
T1_enriched_all = HEC_T1_diff_iso_and_gene_tab[HEC_T1_diff_iso_and_gene_tab$enrich_type=="T1_pre_HSC_enrich",
c("transcript_id","transcript_name","gene_id","gene_name","enrich_type")]
T1_enriched = T1_enriched_all[setdiff(row.names(T1_enriched_all), union(row.names(AEC_enriched),row.names(HEC_enriched))),]

all_sig_iso_and_gene = unique(rbind(AEC_enriched, HEC_enriched, T1_enriched))
AEC_HEC_T1_diff_transcript_exp = transcript_tpm[row.names(all_sig_iso_and_gene),samples]
AEC_HEC_T1_diff_transcript_median_exp = t(apply(AEC_HEC_T1_diff_transcript_exp,1,function(x)
{tapply(x,as.character(sample2stage[names(x),"stage"]),median)}))
all_sig_iso_and_gene_frm = data.frame(all_sig_iso_and_gene, AEC_HEC_T1_diff_transcript_median_exp[row.names(all_sig_iso_and_gene),])


all_sig_gene = unique(all_sig_iso_and_gene[,c("gene_id","gene_name","enrich_type")])
all_sig_gene_id = unique(as.character(all_sig_gene$gene_id))
all_sig_gene_two = as.character(all_sig_gene[duplicated(all_sig_gene$gene_id),"gene_id"])
all_sig_gene_one = setdiff(as.character(all_sig_gene$gene_id),all_sig_gene_two)
all_sig_gene_two_tab = all_sig_gene[all_sig_gene$gene_id %in% all_sig_gene_two,]
all_sig_gene_two_tab = all_sig_gene_two_tab[order(all_sig_gene_two_tab$gene_id),]
all_sig_gene_tpm = gene_tpm[all_sig_gene_id,samples]
all_sig_gene_median_exp = t(apply(all_sig_gene_tpm,1,function(x)
{tapply(x,as.character(sample2stage[names(x),"stage"]),median)}))
all_sig_gene_two_tab_exp = data.frame(all_sig_gene_two_tab, all_sig_gene_median_exp[as.character(all_sig_gene_two_tab$gene_id),])
all_sig_gene_two_tab_exp$max_exp = apply(all_sig_gene_two_tab_exp[,4:6],1,function(x){paste(names(x[order(x,decreasing=T)])[1],"_enrich",sep="")}) 
all_sig_gene_two_tab_exp_filter = all_sig_gene_two_tab_exp[apply(all_sig_gene_two_tab_exp,1,function(x){
x["enrich_type"]==x["max_exp"]}),]
row.names(all_sig_gene_two_tab_exp_filter) = all_sig_gene_two_tab_exp_filter$gene_id
all_sig_gene_two_tab$enrich_type = all_sig_gene_two_tab_exp_filter[as.character(all_sig_gene_two_tab$gene_id),"enrich_type"]
all_sig_gene_two_unique = unique(all_sig_gene_two_tab)
all_sig_gene_tab = rbind(all_sig_gene[all_sig_gene$gene_id %in% all_sig_gene_one,],all_sig_gene_two_unique)
row.names(all_sig_gene_tab) = all_sig_gene_tab$gene_id

all_sig_gene_for_plot = data.frame(all_sig_gene_tab,all_sig_gene_median_exp[as.character(all_sig_gene_tab$gene_id),])

plot_frm = all_sig_gene_for_plot
plot_frm$enrich_type = factor(plot_frm$enrich_type,levels = c("AEC_enrich","HEC_enrich","T1_pre_HSC_enrich"))

library("ggtern")
text_size = 8
line_size = 0.5/1.07

p1<-ggtern(data=plot_frm,aes(x=HEC,y=AEC,z=T1_pre_HSC, color= enrich_type,label = gene_name))+
geom_point(alpha=1,size=2,shape=16)+
#geom_text(color="black",size=2)+
scale_color_manual(values = c("#757575","#689f38","#f44336"))+
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


write.table(plot_frm,"./diff_iso_and_gene/AEC_HEC_T1/AEC_HEC_T1_sig_iso_and_gene_table.txt",sep="\t",quote=F,row.names=F)


candi_gene = c("Runx1","Dnmt3b")
candi_plot_frm = plot_frm[plot_frm$gene_name %in% candi_gene,]
candi_plot_frm$enrich_type = candi_plot_frm$gene_name
candi_plot_frm_all = plot_frm[setdiff(row.names(plot_frm),row.names(candi_plot_frm)),]
candi_plot_frm_all$gene_name = rep("",nrow(candi_plot_frm_all))
candi_plot_frm_all_new = rbind(candi_plot_frm_all,candi_plot_frm)
plot_tab = candi_plot_frm_all_new

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

