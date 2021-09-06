rm(list=ls())
gc()

rmaps_dir = './bed2rMAPS'
######################################################################################
#SE
######################################################################################
dir.create("event_included_exon")

SE = read.table(paste(rmaps_dir,"SE_Events_loc_for_rMAPS.txt",sep='/'),sep="\t",header=T,row.names=1,check.names=F)
SE_event_new = apply(SE,1,function(x){paste(x['chr'],
paste(as.numeric(x['exonStart']),as.numeric(x['exonEnd']),sep='-'),x['strand'],sep=':')})
SE_event_tab = as.data.frame(SE_event_new)

SE_event_frm = data.frame(event_name = row.names(SE_event_tab), 
SE_event_tab, event_type = rep("SE",nrow(SE_event_tab)),check.rows=T)
colnames(SE_event_frm) = c('event_name','alt_exon','event_type')
write.table(SE_event_frm,"./event_included_exon/SE_EventIncludedExon.txt",sep="\t",quote=F,row.names=F)


######################################################################################
#A3SS
######################################################################################

A3SS = read.table(paste(rmaps_dir,"A3SS_Events_loc_for_rMAPS.txt",sep='/'),sep="\t",header=T,row.names=1,check.names=F)
A3SS_pos = A3SS[A3SS[,'strand']=='+',]
A3SS_pos_new = apply(A3SS_pos,1,function(x){paste(x['chr'],
paste(as.numeric(x['longexonStart']),as.numeric(x['shortExonStart'])-1,sep='-'),x['strand'],sep=':')})
A3SS_pos_tab = as.data.frame(A3SS_pos_new)
colnames(A3SS_pos_tab) = 'alt_exon'

A3SS_neg = A3SS[A3SS[,'strand']=='-',]
A3SS_neg_new = apply(A3SS_neg,1,function(x){paste(x['chr'],
paste(as.numeric(x["shortExonEnd"])+1,as.numeric(x["longexonEnd"]),sep="-"),x['strand'],sep=':')})

A3SS_neg_tab = as.data.frame(A3SS_neg_new)
colnames(A3SS_neg_tab) = 'alt_exon'
A3SS_tab = rbind(A3SS_pos_tab, A3SS_neg_tab)

A3SS_tab_frm = data.frame(event_name = row.names(A3SS_tab), A3SS_tab, event_type = rep("A3SS",nrow(A3SS_tab)),check.rows=T)
write.table(A3SS_tab_frm,"./event_included_exon/A3SS_EventIncludedExon.txt",sep="\t",quote=F,row.names=F)

######################################################################################
#A5SS
######################################################################################

A5SS = read.table(paste(rmaps_dir,"A5SS_Events_loc_for_rMAPS.txt",sep='/'),sep="\t",header=T,row.names=1,check.names=F)
A5SS_pos = A5SS[A5SS[,'strand']=='+',]
A5SS_pos_new = apply(A5SS_pos,1,function(x){paste(x[1],
paste(as.numeric(x[6])+1,as.numeric(x[4]),sep="-"),x[2],sep=':')})
A5SS_pos_tab = as.data.frame(A5SS_pos_new)
colnames(A5SS_pos_tab) = 'alt_exon'

A5SS_neg = A5SS[A5SS[,'strand']=='-',]
A5SS_neg_new = apply(A5SS_neg,1,function(x){paste(x[1],
paste(as.numeric(x[3]),as.numeric(x[5])-1,sep="-"),x[2],sep=':')})
A5SS_neg_tab = as.data.frame(A5SS_neg_new)
colnames(A5SS_neg_tab) = 'alt_exon'

A5SS_tab = rbind(A5SS_pos_tab, A5SS_neg_tab)
A5SS_tab_frm = data.frame(event_name = row.names(A5SS_tab), A5SS_tab, event_type = rep('A5SS',nrow(A5SS_tab)),check.rows=T)
write.table(A5SS_tab_frm,"./event_included_exon/A5SS_EventIncludedExon.txt",sep="\t",quote=F,row.names=F)

######################################################################################
#A3SS
######################################################################################
MXE = read.table(paste(rmaps_dir,"MXE_Events_loc_for_rMAPS.txt",sep='/'),sep="\t",header=T,row.names=1,check.names=F)
MXE_alt_exon = apply(MXE,1,function(x){paste(x[1],paste(as.numeric(x[3]),as.numeric(x[4]),sep='-'),x[2],sep=':')})
MXE_tab_frm = data.frame(event_name = names(MXE_alt_exon), alt_exon = MXE_alt_exon, event_type = rep("MXE",length(MXE_alt_exon)))
write.table(MXE_tab_frm,"./event_included_exon/MXE_EventIncludedExon.txt",sep="\t",quote=F,row.names=F)


######################################################################################
#RI
######################################################################################

RI = read.table(paste(rmaps_dir,"RI_Events_loc_for_rMAPS.txt",sep='/'),sep="\t",header=T,row.names=1,check.names=F)
RI_new = apply(RI,1,function(x){paste(x[1],
paste(as.numeric(x[6])+1,as.numeric(x[7])-1,sep="-"),x[2],sep=':')})
RI_tab = as.data.frame(RI_new)
RI_frm = data.frame(event_name = row.names(RI_tab),RI_tab,event_type = rep("RI",nrow(RI_tab)),check.rows=T)
colnames(RI_frm) = c('event_name','alt_exon','event_type')
write.table(RI_frm,"./event_included_exon/RI_EventIncludedExon.txt",sep="\t",quote=F,row.names=F)

######################################################################################
#merge all events together.
######################################################################################

all_event = rbind(SE_event_frm, A3SS_tab_frm, A5SS_tab_frm, MXE_tab_frm, RI_frm)
write.table(all_event,"./event_included_exon/allEventIncludedExon.txt",sep="\t",quote=F,row.names=F)

all_event_bed = t(apply(all_event,1,function(x){unlist(strsplit(x[2],split=":"))}))
all_event_bed2 = t(apply(all_event_bed,1,function(x){unlist(strsplit(x[2],split="-"))}))
score = apply(all_event_bed2,1,function(x){as.numeric(x[2])-as.numeric(x[1])})

all_event_bed_tab = data.frame(chr = all_event_bed[names(score),1],start = all_event_bed2[names(score),1],
end = all_event_bed2[names(score),2],event_name = names(score),score = score,strand = all_event_bed[names(score),3])
write.table(all_event_bed_tab,"./event_included_exon/all_included_region_loc.bed",sep="\t",quote=F,row.names=F)






