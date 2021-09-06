rm(list=ls())
gc()
########################################################################################################
#read event infor and event modality.
########################################################################################################
modality_tab = read.table("../Anchor/EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1,check.names=F)
event_infor = read.table("../merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1,check.names=F)
event_psi = read.table("../merge_PSI/merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,row.names=1,check.names=F)
sample2stage = read.table("../../sample_infor/sample2stage.txt",header=T,sep="\t",stringsAsFactors=F)
row.names(sample2stage) = sample2stage[,1]

event_median_PSI = t(apply(event_psi,1,function(x){tapply(x,sample2stage[names(x),'stage'],function(y){median(as.numeric(y),na.rm=T)})}))

########################################################################################################
#get rMAPS input files.
########################################################################################################
getrMAPS <- function(event, type){
rMAPS_path = "../merge_PSI/bed2rMAPS"
tab = read.table(paste(rMAPS_path, paste(type,"_Events_loc_for_rMAPS.txt",sep=""),sep="/"),sep="\t",row.names=1,header=T)
intersect_event = intersect(row.names(tab), event)
out_tab = tab[intersect_event,]
return(out_tab)
}

##################################################################################
#get T1 signature events.
##################################################################################
getSigEvents <- function(mod_tab, stage, mod_type = "included", eventInfor = event_infor, eventPSI, outPath){
pos_stage = stage
neg_stage = setdiff(colnames(mod_tab),pos_stage)
#------filter NULL line in modality data frame.
median_PSI = eventPSI
median_PSI[is.na(median_PSI)]=0
colnames(median_PSI) = paste(colnames(eventPSI),"_median_PSI",sep="")
non_null_tab = mod_tab[apply(mod_tab,1,function(x){length(which(x==""))!=length(x)}),c(pos_stage, neg_stage)]
non_null_frm = data.frame(non_null_tab,eventInfor[row.names(non_null_tab),],
median_PSI[row.names(non_null_tab),],check.names=F)
non_null_frm$delta_psi = apply(non_null_frm,1,function(x){
as.numeric(x[paste(pos_stage,"_median_PSI",sep="")])-
median(as.numeric(x[paste(neg_stage,"_median_PSI",sep="")],na.rm=T))})

#------get upregulated events.
pos_sig_tab = non_null_frm[non_null_frm[,pos_stage]==mod_type,]
pos_sig_event = pos_sig_tab[pos_sig_tab[,neg_stage]!=mod_type,]
pos_sig_event_filter = pos_sig_event

#------get downregulated events.
neg_sig_tab = non_null_frm[non_null_frm[,pos_stage] != mod_type,]
neg_sig_event = neg_sig_tab[neg_sig_tab[,neg_stage]==mod_type,]
neg_sig_event_filter = neg_sig_event

#------get background events.
back_sig_tab = non_null_frm[apply(non_null_frm,1,function(x){x[pos_stage] == x[neg_stage]}),]
back_sig_event = back_sig_tab[back_sig_tab[,pos_stage]!="multimodal",]
back_sig_event_filter = back_sig_event[abs(back_sig_event$delta_psi)<0.05,]

#------merge 3 set together.
three_set_modality_tab = paste(outPath,"three_set_modality_tab",sep='/')
dir.create(three_set_modality_tab, recursive = T, showWarnings=F)
write.table(pos_sig_event_filter, file.path(three_set_modality_tab,paste(stage, mod_type, "events_postive_set.txt",sep="_")),sep="\t",quote=F)
write.table(neg_sig_event_filter, file.path(three_set_modality_tab,paste(stage, mod_type, "events_negative_set.txt",sep="_")),sep="\t",quote=F)
write.table(back_sig_event_filter, file.path(three_set_modality_tab,paste(stage, mod_type, "events_background_set.txt",sep="_")),sep="\t",quote=F)

#------get 3 set for rMAPS.

event_type = c("SE","A3SS","A5SS","MXE","RI")
for(i in event_type){
outdir = paste(outPath,'signature_events_for_rMAPS',i,sep = '/')
dir.create(outdir, recursive = T, showWarnings=F)
#------postive set
up_rmaps_event = row.names(pos_sig_event_filter[pos_sig_event_filter$event_type == i,])
up_rmaps_tab = getrMAPS(event = up_rmaps_event, type = i)
write.table(up_rmaps_tab,file.path(outdir,paste(stage, mod_type, "events_postive_set.txt",sep="_")),sep="\t",quote=F,row.names=F,col.names=F)
#------negative set
down_rmaps_event = row.names(neg_sig_event_filter[neg_sig_event_filter$event_type == i,])
down_rmaps_tab = getrMAPS(event = down_rmaps_event, type = i)
write.table(down_rmaps_tab,file.path(outdir,paste(stage, mod_type, "events_negative_set.txt",sep="_")),sep="\t",quote=F,row.names=F,col.names=F)
#------postive set
back_rmaps_event = row.names(back_sig_event_filter[back_sig_event_filter$event_type == i,])
back_rmaps_tab = getrMAPS(event = back_rmaps_event, type = i)
write.table(back_rmaps_tab,file.path(outdir,paste(stage, mod_type, "events_background_set.txt",sep="_")),sep="\t",quote=F,row.names=F,col.names=F)
}
}

##################################################################################

#------AEC HEC
AEC_HEC_sig = modality_tab[,c("AEC","HEC")]
AEC_HEC_sig_outDir = "./rMAPS_enrichment/AEC_HEC"
AEC_HEC_median_PSI = event_median_PSI[,c("AEC","HEC")]
AEC_HEC_sig_in_rMAPS = getSigEvents(mod_tab = AEC_HEC_sig, stage="HEC", eventPSI = AEC_HEC_median_PSI, outPath = AEC_HEC_sig_outDir)

#------HEC T1 pre-HSC
HEC_T1_sig = modality_tab[,c("HEC","T1_pre_HSC")]
HEC_T1_sig_outDir = "./rMAPS_enrichment/HEC_T1"
HEC_T1_median_PSI = event_median_PSI[,c("HEC","T1_pre_HSC")]
HEC_T1_sig_in_rMAPS = getSigEvents(mod_tab = HEC_T1_sig, stage="T1_pre_HSC", eventPSI = HEC_T1_median_PSI, outPath = HEC_T1_sig_outDir)

#------AEC T1 pre-HSC
AEC_T1_sig = modality_tab[,c("AEC","T1_pre_HSC")]
AEC_T1_sig_outDir = "./rMAPS_enrichment/AEC_T1"
AEC_T1_median_PSI = event_median_PSI[,c("AEC","T1_pre_HSC")]
AEC_T1_sig_in_rMAPS = getSigEvents(mod_tab = AEC_T1_sig, stage="T1_pre_HSC", eventPSI = AEC_T1_median_PSI, outPath = AEC_T1_sig_outDir)




