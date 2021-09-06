rm(list=ls())
gc()

HEC_T1_in = read.table("../Anchor/AEC_HEC_T1/HEC_T1_included__AEC_not_events.txt",sep="\t",header=T,row.names=1,check.names=F)
T1_only_in = read.table("../Anchor/AEC_HEC_T1/T1_included_AEC_HEC_not_included_events.txt",sep="\t",header=T,row.names=1,check.names=F)

RBPmap_loc = read.table("../merge_PSI/RBPMap_loc/with_300bp_intron/SE_Events_for_RBPmap.txt",sep="\t",header=T,row.names=1,check.names=F)

dir.create("RBPmap_res/Srsf2_loc",recursive = T)
HEC_T1_in_SE = HEC_T1_in[HEC_T1_in$event_type=="SE",]
HEC_T1_in_SE_RBPmap = as.character(RBPmap_loc[intersect(row.names(HEC_T1_in_SE),row.names(RBPmap_loc)),])
write.table(HEC_T1_in_SE_RBPmap,"./RBPmap_res/Srsf2_loc/HEC_T1_in_SE_loc.txt",sep="\t",row.names=F,col.names=F,quote=F)


T1_only_in_SE = T1_only_in[T1_only_in$event_type=="SE",]
T1_only_in_SE_RBPmap = as.character(RBPmap_loc[intersect(row.names(T1_only_in_SE),row.names(RBPmap_loc)),])
write.table(T1_only_in_SE_RBPmap,"./RBPmap_res/Srsf2_loc/T1_only_in_SE_loc.txt",sep="\t",row.names=F,col.names=F,quote=F)

