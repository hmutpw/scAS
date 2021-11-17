#---
sample2stage <- read_delim("D:/YJ/YJ-AS/NewAnalysis_AGM/Nanopore_seq/sample_infor/sample2stage.txt",delim = "\t")
transcript_infor <- read_delim(file = "D:/YJ/YJ-AS/NewAnalysis_AGM/ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",delim = "\t")
gene_infor <- read_delim(file = "D:/YJ/YJ-AS/NewAnalysis_AGM/ref_genome/gencode.vM22.gene.infor.rm.version.tsv",delim = "\t")

transcript_mat <- read_delim("D:/YJ/YJ-AS/NewAnalysis_AGM/Nanopore_seq/expression_data/merged_exp/normalized_transcript_RPT10K.tsv",delim ="\t", col_names =T)
gene_mat <- read_delim("D:/YJ/YJ-AS/NewAnalysis_AGM/Nanopore_seq/expression_data/merged_exp/normalized_gene_RPG10K.tsv",delim ="\t", col_names =T)

######merge data.frame for plot
merge_frm <- function(gene, iso_exp = transcript_mat, gene_exp = gene_mat ,iso_infor = transcript_infor, g_infor=gene_infor, sample2Stage=sample2stage){
  isoform_infor <- iso_infor %>% filter(gene_name == gene)
  gene_infor <- g_infor %>% filter(gene_name == gene)
  candidate_gene_exp <- gene_exp %>% filter(gene_id %in% gene_infor$gene_id) %>% gather(-gene_id,key="sample",value="counts") %>% left_join(.,y=gene_infor,by=c("gene_id"="gene_id")) %>% left_join(.,y=sample2Stage, by=c("sample"="sample")) %>% select(gene_id,sample,counts,gene_name,stage)
  candidate_isoform_exp <- iso_exp %>% filter(transcript_id %in% isoform_infor$transcript_id) %>% gather(-transcript_id,key="sample",value="counts")%>% left_join(.,y=isoform_infor[,c(1,2,5)],by=c("transcript_id"="transcript_id")) %>% left_join(.,y=sample2Stage, by=c("sample"="sample")) %>% select(transcript_id,sample,counts,transcript_name,stage)
colnames(candidate_isoform_exp) <- c("gene_id","sample","counts","gene_name","stage")
  rbind(candidate_gene_exp, candidate_isoform_exp)
}
#---Runx1
Runx1_exp <- merge_frm(gene = "Runx1")
#---Polr2m
Polr2m_exp <-  merge_frm(gene = "Polr2m")
#---Sf3b1
Sf3b1_exp <- merge_frm(gene = "Sf3b1")
#---Dnmt3b
Dnmt3b_exp <- merge_frm(gene = "Dnmt3b")
#---Tulp4
Tulp4_exp <- merge_frm(gene = "Tulp4")
#---E130309D02Rik
E130309D02Rik_exp <- merge_frm(gene = "E130309D02Rik")
#---Pisd
Pisd_exp <- merge_frm(gene = "Pisd")
#---Sec31a
Sec31a_exp <- merge_frm(gene = "Sec31a")
#---Tmpo
Tmpo_exp <-  merge_frm(gene = "Tmpo")
#---Ntmt1
Ntmt1_exp <- merge_frm(gene = "Ntmt1")
#---Clk1
Clk1_exp <- merge_frm(gene = "Clk1")


text_size = 8
line_size = 0.5/1.07
library(ggplot2)
#------set plot theme
mytheme = theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA,color=NA),
strip.background = element_rect(fill=NA,color=NA),
strip.text = element_text(size=text_size,color="black"),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
panel.grid.minor = element_line(size=NA, color = NA),
axis.ticks = element_line(size=line_size,color="black"),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

plot_frm <- Dnmt3b_exp
plot_frm <- plot_frm %>% filter(stage %in% c("AEC","HEC")) %>% filter(gene_name %in% c("Dnmt3b","Dnmt3b-201","Dnmt3b-206"))


p <- ggplot(plot_frm, mapping = aes(x=stage,y=counts, color=stage))+
#geom_violin(position = "dodge", na.rm = TRUE, scale = "width",width = 0.9)+
geom_boxplot(width=0.25,fill=NA, size=0.2,outlier.colour=NA, alpha = 0.8)+
geom_jitter(width=0.2,size=0.8,shape=16,fill="white", alpha = 0.8)+
  scale_color_manual(values = c("#1e88e5","#afb42b"))+
labs(x="",y='RPT10K')+
mytheme+
facet_wrap(.~gene_name,nrow=1)

p

