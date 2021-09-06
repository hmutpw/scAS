mergeModality <- function(stage,inpath = inDir){
modality = read.table(paste(inpath,paste(stage,"_all_events_modality.txt",sep=''),sep="/"),sep="\t",header=T)
colnames(modality) = c("event_name",stage)
return(modality)
}
inDir = "./EventModality/Events2Modality"

stages = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
mergedModality = mergeModality(stage = stages[1])

for(i in stages[-1]){
tab = mergeModality(stage = i)
mergedModality = merge(mergedModality,tab,by = "event_name",all=T)
}

outPath = "./EventModality"
write.table(mergedModality,file.path(outPath,"7_stage_event_modality.txt"),sep="\t",quote=F,row.names=F,na="")


#------add gene infor to modality.
modality_tab = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T)
row.names(modality_tab) = modality_tab$event_name
event_infor = read.table("../merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1)

modality_tab_out = data.frame(modality_tab,event_infor[row.names(modality_tab),],check.names = F)
write.table(modality_tab_out,"./EventModality/7_stage_event_modality_with_gene_infor.txt",sep="\t",quote=F,row.names=F,na="")







