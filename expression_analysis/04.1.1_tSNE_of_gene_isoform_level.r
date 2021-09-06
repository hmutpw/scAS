rm(list=ls())
gc()
#install.packages("Seurat")
library(Seurat)


transcript_infor = read.table("../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
sample2stage = read.table("../sample_infor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage$sample
stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
samples = as.character(sample2stage[sample2stage$stage %in% stage,"sample"])

transcript_tpm = read.table("../expression_data/merged_exp/7_stage_transcript_tpm.txt",header=T,row.names=1,check.names=F)
gene_tpm = read.table("../expression_data/merged_exp/7_stage_gene_tpm.txt",header=T,row.names=1,check.names=F)
#ERCC_spike_in = read.table("../expression_data/merged_exp/gene_ERCC_counts.txt",sep="\t",header=T,row.names=1,check.names=F)

counts_mat = transcript_tpm[,samples]
#counts_mat = gene_tpm[,samples]

counts_mat_filter = counts_mat[apply(counts_mat,1,function(x){length(x[which(x>1)])>=2}),]
#------Step1. prepare data for analysis
pbmc.data = counts_mat_filter
pbmc <- CreateSeuratObject(counts = pbmc.data)
pbmc

orig.ident = factor(as.character(sample2stage[samples,'stage']),levels = stage)
pbmc@meta.data[,"orig.ident"]=orig.ident
pbmc@active.ident = orig.ident
names(pbmc@active.ident) = row.names(pbmc@meta.data)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

#------Step2. Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#------Step3. Detection of variable genes across the single cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)

#------Step4. Scaling the data and removing unwanted sources of variation.
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

#------Step5. Perform linear dimensional reduction
pbmc <- RunPCA(pbmc,npcs=100)
DimPlot(pbmc, reduction = "pca")

pbmc <- JackStraw(pbmc, num.replicate = 10, dims = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:100)
JackStrawPlot(pbmc, dims = 1:100)
ElbowPlot(pbmc, ndims = 100)

colors = c("#2080DF","#80DF20","#DF2020","#20DF80","#DF8020","#8020DF","#DF2080")

pbmc <- RunTSNE(object = pbmc, dims.use = 1:30)
DimPlot(pbmc, reduction = "tsne",label = TRUE)

#reticulate::py_install(packages = 'umap-learn')

pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc, reduction = "umap",label = TRUE,cols = colors, pt.size =1,label.size = 4)

umap_res = pbmc@reductions$umap@cell.embeddings
write.table(umap_res,"iso_level_umap_res.txt",sep="\t",quote=F)

dir.create("seurat_output")
saveRDS(pbmc, file = "./seurat_output/isoform_level.rds")






