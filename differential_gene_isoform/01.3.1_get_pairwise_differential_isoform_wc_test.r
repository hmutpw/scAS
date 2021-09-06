rm(list=ls())
gc()
####################################################################################
#------get transcript & sample infor.
####################################################################################
library(data.table)
transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",row.names=1,header=T,stringsAsFactors=F)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage)= sample2stage$sample
stages = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = row.names(sample2stage[sample2stage$stage %in% stages,])

#------get gene TPM.
exp_dir = "../expression_data/merged_exp"
transcript_tpm = fread(file.path(exp_dir, "7_stage_transcript_TPM.txt"),sep="\t",header=T,data.table=F)
tpm_mat = transcript_tpm[,-1]
row.names(tpm_mat) = transcript_tpm$transcript_id


########################################################################################################
#------wilcox.test
########################################################################################################
deTest <- function(case, control, min.cell = 5, min.exp = 1){
#------get case and control sample name.
case_samp = colnames(case)
control_samp = colnames(control)
#------filter gene with more than min.cell.
case_exp = case[apply(case,1,function(x){length(which(x>min.exp))>=min.cell}),]
control_exp = control[apply(control,1,function(x){length(which(x>min.exp))>=min.cell}),]
#------merge exp matrix with row names.
exp_mat = merge(case_exp, control_exp, by = 0, all=T, sort = F)
exp_mat_filter = exp_mat[,-1]
row.names(exp_mat_filter) = exp_mat[,1]
exp_mat_filter[is.na(exp_mat_filter)] = 0
#------perform wilcoxon rank sum test.
p.value = apply(exp_mat_filter,1,function(x){wilcox.test(as.numeric(x[case_samp]),as.numeric(x[control_samp]),exact=F)$p.value})
FDR = p.adjust(p.value, method = "BH")

median_Set1 = apply(exp_mat_filter[,case_samp],1,median)
median_Set2 = apply(exp_mat_filter[,control_samp],1,median)
logFC = log2((median_Set2+1e-2)/(median_Set1+1e-2))

outFrame = data.frame(median_Set1 = median_Set1[names(median_Set1)], 
median_Set2 = median_Set2[names(median_Set1)], 
logFC = logFC[names(median_Set1)], 
p_value = p.value[names(median_Set1)], 
FDR = FDR[names(median_Set1)], check.rows = T)
outFrameSort = outFrame[order(outFrame$FDR),]
return(outFrameSort)
}

########################################################################################################
#perform Test
########################################################################################################

getSigGenes <- function(stage1, stage2, outpath = outPath, cutoff = 0.05, gene_tpm = tpm_mat, sampleStage = sample2stage, GeneInfor = transcript_infor){
sample1 = as.character(sampleStage[sampleStage$stage %in% stage1,"sample"])
sample2 = as.character(sampleStage[sampleStage$stage %in% stage2,"sample"])

Set1 = gene_tpm[,sample1]
Set2 = gene_tpm[,sample2]
DE_gene = deTest(case = Set1, control = Set2, min.cell = 5)
intersect_gene = intersect(row.names(DE_gene), row.names(GeneInfor))
DE_gene_frm = data.frame(transcript_id = intersect_gene,
GeneInfor[intersect_gene,],
DE_gene[intersect_gene,], check.rows=T,check.names=F)
colnames(DE_gene_frm)[7:8] = paste(c(stage1,stage2),'median_TPM',sep='_')
all_res_path = paste(outpath,"all_res",sep='/')
dir.create(all_res_path, recursive = TRUE, showWarnings = F)
write.table(DE_gene_frm, paste(all_res_path,paste(stage1,stage2,'test_all_res.txt',sep='_'),sep='/'),sep = "\t",row.names=F, quote=F)

sig_gene = DE_gene_frm[DE_gene_frm$FDR < cutoff,]
sig_gene_lfc = sig_gene[abs(sig_gene$logFC)>=1,]
sig_gene_high = sig_gene_lfc[apply(sig_gene_lfc[,grep("_median_TPM$",colnames(sig_gene_lfc))],1,function(x){length(x[x>1])>0}),]
diff_res_path = paste(outpath,"diff_res",sep='/')
dir.create(diff_res_path, recursive = TRUE, showWarnings = F)
write.table(sig_gene_high,paste(diff_res_path,paste(stage1,stage2,'sig_diff_res.txt',sep='_'),sep='/'),sep = "\t",row.names=F, quote=F)
return(sig_gene_lfc)
}


outPath = "./diff_res/transcript_level/pairwise_test"
AEC_HEC_sig = getSigGenes(stage1 = "AEC", stage2 = "HEC")
HEC_T1_pre_HSC_sig = getSigGenes(stage1 = "HEC", stage2 = "T1_pre_HSC")
AEC_T1_pre_HSC_sig = getSigGenes(stage1 = "AEC", stage2 = "T1_pre_HSC")
T1_pre_HSC_T2_pre_HSC_sig = getSigGenes(stage1 = "T1_pre_HSC", stage2 = "T2_pre_HSC")
T2_pre_HSC_E12_sig = getSigGenes(stage1 = "T2_pre_HSC", stage2 = "E12")
E12_E14_sig = getSigGenes(stage1 = "E12", stage2 = "E14")
E14_Adult_HSC_sig = getSigGenes(stage1 = "E14", stage2 = "Adult_HSC")




