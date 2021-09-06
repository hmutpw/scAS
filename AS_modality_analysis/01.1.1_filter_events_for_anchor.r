PSI_path = "../merge_PSI/merged_PSI/PSI_each_stage"
out_Path  = "./PSIforAnchor"
dir.create(out_Path)
filterPSI <- function(Stage, PSINumCutoff = 5, inPath = PSI_path, outPath = out_Path){
psi = read.table(file.path(inPath,paste(Stage,'_PSI_Tab.txt',sep='')),sep="\t",header=T,row.names=1)
event = psi[apply(psi,1,function(x){length(which(!is.na(x)))>=PSINumCutoff}),]
event_t = t(event)
write.table(event_t,paste(outPath,paste(Stage,"_PSI_Filter_",PSINumCutoff,"_Tab.txt",sep=""),sep="/"),sep="\t",quote=F)
}

stages = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
for(i in stages){
filterPSI(Stage = i)
}