rm(list=ls())
sample2stage = read.table("../../sample_infor/sample2stage.txt",sep="\t",header=T,stringsAsFactors=F)
row.names(sample2stage) = sample2stage$sample

psi_dir = "../miso_out_res/summary_output"
list.files(psi_dir)
library(data.table)

##########################
#get PSI for each sample.
##########################
getPSI <- function(sample, inputPath){
SE = fread(file.path(inputPath,sample,'SE/summary/SE.miso_summary'),header=T,sep="\t",data.table=F)
A3SS = fread(file.path(inputPath,sample,'A3SS/summary/A3SS.miso_summary'),header=T,sep="\t",data.table=F)
A5SS = fread(file.path(inputPath,sample,'A5SS/summary/A5SS.miso_summary'),header=T,sep="\t",data.table=F)
MXE = fread(file.path(inputPath,sample,'MXE/summary/MXE.miso_summary'),header=T,sep="\t",data.table=F)
RI = fread(file.path(inputPath,sample,'RI/summary/RI.miso_summary'),header=T,sep="\t",data.table=F)
all_event = rbind(SE, A3SS, A5SS, MXE, RI)
psi_tab = all_event[,c("event_name","miso_posterior_mean")]
colnames(psi_tab) = c("event_name",sample)
return(psi_tab)
}
##########################
#merge PSI for each stage.
##########################
getStagePSI <- function(stage, input_dir = inPath, output_dir = outPath, sampleStage = sample2stage){
samplelist = as.character(sampleStage[sampleStage$stage==stage,"sample"])
#------get PSI tab for each stage.
PSITab = getPSI(sample = samplelist[1], inputPath = input_dir)
for(i in samplelist[-1]){
psitab = getPSI(sample = i, inputPath = input_dir)
PSITab = merge(x = PSITab, y = psitab, by = "event_name", all=TRUE)
}
out_stage_dir = file.path(output_dir,"PSI_each_stage")
dir.create(out_stage_dir, recursive = T)
fwrite(PSITab[order(PSITab[,"event_name"]),], file.path(out_stage_dir, paste(stage,'_PSI_Tab.txt',sep='')),sep="\t",quote=F,row.names=F)
return(PSITab)
}

inPath = psi_dir
outPath = "./merged_PSI"

all_stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
#running merge PSI steps, it takes few minutes.
all_PSI = getStagePSI(stage = all_stage[1])
for(i in all_stage[-1]){
psi = getStagePSI(stage = i)
all_PSI = merge(x = all_PSI, y = psi, by = "event_name", all=TRUE)
}

fwrite(all_PSI, file.path(outPath,'all_PSI_Tab.txt'),sep="\t",quote=F,row.names=F)

all_PSI_frm = fread(file.path(outPath,'all_PSI_Tab.txt'),sep="\t",header=T,data.table=F)
all_PSI = all_PSI_frm[,-1]
all_PSI_num = apply(all_PSI,2,function(x){length(x[which(!is.na(x))])})
write.table(all_PSI_num, 'all_sample_event_number.txt',sep="\t",quote=F)


