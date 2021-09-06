rmaps_path = "./bed2rMAPS"
rbpmap_path = "./RBPMap_loc/with_300bp_intron"
dir.create(rbpmap_path,recursive = T)
############
#SE
############
SE_loc = read.table(paste(rmaps_path,"SE_Events_loc_for_rMAPS.txt",sep="/"),sep="\t",header=T,row.names=1)
SE_loc$up_intron_len = apply(SE_loc,1,function(x){
as.numeric(x["exonStart"])-as.numeric(x["firstExonEnd"])})
SE_loc$down_intron_len = apply(SE_loc,1,function(x){
as.numeric(x["secondExonStart"])-as.numeric(x["exonEnd"])})

intron_length = 300

SE_loc_for_rbpmap = apply(SE_loc,1,function(x){paste(x["chr"],
paste(as.numeric(x["exonStart"])-intron_length,
as.numeric(x["exonEnd"])+intron_length
,sep='-'),x["strand"],sep=':')})
SE_loc_for_rbpmap_tab = data.frame(event_name = names(SE_loc_for_rbpmap),
loc = SE_loc_for_rbpmap)

write.table(SE_loc_for_rbpmap_tab,paste(rbpmap_path, "SE_Events_for_RBPmap.txt",sep="/"),sep="\t",
quote=F,row.names=F)

############
#A3SS
############
A3SS_loc = read.table(paste(rmaps_path,"A3SS_Events_loc_for_rMAPS.txt",sep="/"),header=T,row.names=1,sep="\t")
A3SS_loc_pos = A3SS_loc[A3SS_loc$strand=="+",]
A3SS_loc_neg = A3SS_loc[A3SS_loc$strand=="-",]
A3SS_loc_pos_for_rbpmap = apply(A3SS_loc_pos,1,function(x){paste(x["chr"],
paste(as.numeric(x["longexonStart"])-intron_length,as.numeric(x["longexonEnd"])
,sep='-'),x["strand"],sep=':')})
A3SS_loc_neg_for_rbpmap = apply(A3SS_loc_neg,1,function(x){paste(x["chr"],
paste(as.numeric(x["longexonStart"]),as.numeric(x["longexonEnd"])+intron_length
,sep='-'),x["strand"],sep=':')})
A3SS_loc_for_rbpmap = c(A3SS_loc_pos_for_rbpmap, A3SS_loc_neg_for_rbpmap)
A3SS_loc_for_rbpmap_tab = data.frame(event_name = names(A3SS_loc_for_rbpmap),
loc = A3SS_loc_for_rbpmap)

write.table(A3SS_loc_for_rbpmap_tab,paste(rbpmap_path, "A3SS_Events_for_RBPmap.txt",sep="/"),sep="\t",
quote=F,row.names=F)

############
#A5SS
############
A5SS_loc = read.table(paste(rmaps_path,"A5SS_Events_loc_for_rMAPS.txt",sep="/"),header=T,row.names=1,sep="\t")
A5SS_loc_pos = A5SS_loc[A5SS_loc$strand=="+",]
A5SS_loc_neg = A5SS_loc[A5SS_loc$strand=="-",]
A5SS_loc_pos_for_rbpmap = apply(A5SS_loc_pos,1,function(x){paste(x["chr"],
paste(as.numeric(x["longexonStart"]),as.numeric(x["longexonEnd"])+intron_length
,sep='-'),x["strand"],sep=':')})
A5SS_loc_neg_for_rbpmap = apply(A5SS_loc_neg,1,function(x){paste(x["chr"],
paste(as.numeric(x["longexonStart"])-intron_length,as.numeric(x["longexonEnd"])
,sep='-'),x["strand"],sep=':')})
A5SS_loc_for_rbpmap = c(A5SS_loc_pos_for_rbpmap, A5SS_loc_neg_for_rbpmap)

A5SS_loc_for_rbpmap_tab = data.frame(event_name = names(A5SS_loc_for_rbpmap),
loc = A5SS_loc_for_rbpmap)

write.table(A5SS_loc_for_rbpmap_tab,paste(rbpmap_path, "A5SS_Events_for_RBPmap.txt",sep="/"),sep="\t",
quote=F,row.names=F)


all_events_loc = rbind(SE_loc_for_rbpmap_tab, A3SS_loc_for_rbpmap_tab, A5SS_loc_for_rbpmap_tab)



