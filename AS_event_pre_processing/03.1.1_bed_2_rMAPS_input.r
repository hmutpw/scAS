bedpath = "./event_bed"
list.files(bedpath)
########################################################################################
#SE
########################################################################################
library(data.table)
#SE formate for rMAPS 2.0.
#chr	strand	exonStart	exonEnd	firstExonStart	firstExonEnd	secondExonStart	secondExonEnd
#
#read SE bed files from step1.
SE_bed = fread(file.path(bedpath,"SE.vM22.exon.bed"),sep="\t",header=F,data.table=F,stringsAsFactors=F)
SE_upstream_bed = SE_bed[SE_bed[,5]=="up",]
row.names(SE_upstream_bed) = SE_upstream_bed[,4]
SE_alternative_bed = SE_bed[SE_bed[,5]=="alt",]
row.names(SE_alternative_bed) = SE_alternative_bed[,4]
SE_downstream_bed = SE_bed[SE_bed[,5]=="down",]
row.names(SE_downstream_bed) = SE_downstream_bed[,4]

SE_event = row.names(SE_alternative_bed)
#events for rMAPS input.
SE_Tab_rMAPS = data.frame(SE_alternative_bed[SE_event,c(1,6,2,3)],
SE_upstream_bed[SE_event,2:3],SE_downstream_bed[SE_event,2:3])
colnames(SE_Tab_rMAPS) = c('chr','strand','exonStart','exonEnd','firstExonStart',
'firstExonEnd','secondExonStart','secondExonEnd')

SE_Tab_for_rMAPS = SE_Tab_rMAPS[order(row.names(SE_Tab_rMAPS)),]
dir.create("./bed2rMAPS")
fwrite(SE_Tab_for_rMAPS,"./bed2rMAPS/SE_Events_loc_for_rMAPS.txt",sep="\t",row.names=T,quote=F)

########################################################################################
#MXE
########################################################################################
#MXE formate for rMAPS 2.0.
#chr	strand	1stExonStart	1stExonEnd	2ndExonStart	2ndExonEnd	upstreamExonStart	upstreamExonEnd	downstreamExonStart	downstreamExonEnd
#
#read MXE bed files from step1.
MXE_bed = fread(file.path(bedpath,"MXE.vM22.exon.bed"),sep="\t",header=F,data.table=F,stringsAsFactors=F)
MXE_upstream_bed = MXE_bed[MXE_bed[,5]=="up",]
row.names(MXE_upstream_bed) = MXE_upstream_bed[,4]
MXE_alternative1_bed = MXE_bed[MXE_bed[,5]=="alt1",]
row.names(MXE_alternative1_bed) = MXE_alternative1_bed[,4]
MXE_alternative2_bed = MXE_bed[MXE_bed[,5]=="alt2",]
row.names(MXE_alternative2_bed) = MXE_alternative2_bed[,4]
MXE_downstream_bed = MXE_bed[MXE_bed[,5]=="down",]
row.names(MXE_downstream_bed) = MXE_downstream_bed[,4]

# '+' strand.
MXE_pos_events = row.names(MXE_alternative1_bed[MXE_alternative1_bed[,6]=="+",])
MXE_pos_Tab = data.frame(MXE_alternative1_bed[MXE_pos_events,c(1,6,2,3)],
MXE_alternative2_bed[MXE_pos_events,c(2,3)],
MXE_upstream_bed[MXE_pos_events,2:3],
MXE_downstream_bed[MXE_pos_events,2:3])
colnames(MXE_pos_Tab) = c('chr','strand','1stExonStart','1stExonEnd','2ndExonStart','2ndExonEnd',
'upstreamExonStart','upstreamExonEnd','downstreamExonStart','downstreamExonEnd')

# '-' strand.
MXE_neg_events = row.names(MXE_alternative1_bed[MXE_alternative1_bed[,6]=="-",])
MXE_neg_Tab = data.frame(MXE_alternative2_bed[MXE_neg_events,c(1,6,2,3)],
MXE_alternative1_bed[MXE_neg_events,c(2,3)],
MXE_upstream_bed[MXE_neg_events,2:3],
MXE_downstream_bed[MXE_neg_events,2:3])
colnames(MXE_neg_Tab) = c('chr','strand','1stExonStart','1stExonEnd','2ndExonStart','2ndExonEnd',
'upstreamExonStart','upstreamExonEnd','downstreamExonStart','downstreamExonEnd')

MXE_Tab_rMAPS = rbind(MXE_pos_Tab, MXE_neg_Tab)
MXE_Tab_for_rMAPS = MXE_Tab_rMAPS[order(row.names(MXE_Tab_rMAPS)),]
write.table(MXE_Tab,"./bed2rMAPS/MXE_Events_loc_for_rMAPS.txt",sep="\t",quote=F)

########################################################################################
#A5SS
########################################################################################
#A5SS formate for rMAPS 2.0.
#chr	strand	longexonStart	longexonEnd	shortExonStart	shortExonEnd	flankingExonStart	flankingExonEnd
#
#read A5SS bed files from step1.
A5SS_bed = fread(file.path(bedpath,"A5SS.vM22.exon.bed"),sep="\t",header=F,data.table=F,stringsAsFactors=F)
A5SS_included_bed = A5SS_bed[A5SS_bed[,5]=="alt_in",]
row.names(A5SS_included_bed) = A5SS_included_bed[,4]
A5SS_excluded_bed = A5SS_bed[A5SS_bed[,5]=="alt_out",]
row.names(A5SS_excluded_bed) = A5SS_excluded_bed[,4]
A5SS_downstream_bed = A5SS_bed[A5SS_bed[,5]=="down",]
row.names(A5SS_downstream_bed) = A5SS_downstream_bed[,4]

A5SS_Tab = data.frame(A5SS_included_bed[,c(1,6,2,3)],
A5SS_excluded_bed[,c(2,3)],A5SS_downstream_bed[,c(2,3)],check.rows=T)
colnames(A5SS_Tab) = c('chr','strand','longexonStart','longexonEnd','shortExonStart','shortExonEnd','flankingExonStart','flankingExonEnd')
A5SS_Tab_for_rMAPS = A5SS_Tab[order(row.names(A5SS_Tab)),]

write.table(A5SS_Tab_for_rMAPS,"./bed2rMAPS/A5SS_Events_loc_for_rMAPS.txt",sep="\t",quote=F)

########################################################################################
#A3SS
########################################################################################
#A3SS formate for rMAPS 2.0.
#chr	strand	longexonStart	longexonEnd	shortExonStart	shortExonEnd	flankingExonStart	flankingExonEnd
#
#read A3SS bed files from step1.
A3SS_bed = fread(file.path(bedpath,"A3SS.vM22.exon.bed"),sep="\t",header=F,data.table=F,stringsAsFactors=F)
A3SS_included_bed = A3SS_bed[A3SS_bed[,5]=="alt_in",]
row.names(A3SS_included_bed) = A3SS_included_bed[,4]
A3SS_excluded_bed = A3SS_bed[A3SS_bed[,5]=="alt_out",]
row.names(A3SS_excluded_bed) = A3SS_excluded_bed[,4]
A3SS_upstream_bed = A3SS_bed[A3SS_bed[,5]=="up",]
row.names(A3SS_upstream_bed) = A3SS_upstream_bed[,4]

A3SS_Tab = data.frame(A3SS_included_bed[,c(1,6,2,3)],A3SS_excluded_bed[,c(2,3)],
A3SS_upstream_bed[,c(2,3)],check.rows=T)
colnames(A3SS_Tab) = c('chr','strand','longexonStart','longexonEnd','shortExonStart',
'shortExonEnd','flankingExonStart','flankingExonEnd')
A3SS_Tab_for_rMAPS = A3SS_Tab[order(row.names(A3SS_Tab)),]

write.table(A3SS_Tab_for_rMAPS,"./bed2rMAPS/A3SS_Events_loc_for_rMAPS.txt",sep="\t",quote=F)

########################################################################################
#RI
########################################################################################
#RI formate for rMAPS 2.0.
#chr	strand	riExonStart	riExonEnd	upstreamExonStart	upstreamExonEnd	downstreamExonStart	downstreamExonEnd
#
#read RI bed files from step1.
RI_bed = fread(file.path(bedpath,"RI.vM22.exon.bed"),sep="\t",header=F,data.table=F,stringsAsFactors=F)
RI_upstream_bed = RI_bed[RI_bed[,5]=="up",]
row.names(RI_upstream_bed) = RI_upstream_bed[,4]
RI_alternative_bed = RI_bed[RI_bed[,5]=="alt",]
row.names(RI_alternative_bed) = RI_alternative_bed[,4]
RI_downstream_bed = RI_bed[RI_bed[,5]=="down",]
row.names(RI_downstream_bed) = RI_downstream_bed[,4]

#'+' strand events for rMAPS input.
RI_Tab_rMAPS = data.frame(RI_alternative_bed[,c(1,6,2,3)],
RI_upstream_bed[row.names(RI_alternative_bed),2:3],
RI_downstream_bed[row.names(RI_alternative_bed),2:3])
colnames(RI_Tab_rMAPS) = c('chr','strand','riExonStart','riExonEnd','upstreamExonStart',
'upstreamExonEnd','downstreamExonStart','downstreamExonEnd')
RI_Tab_for_rMAPS = RI_Tab_rMAPS[order(row.names(RI_Tab_rMAPS)),]
write.table(RI_Tab_for_rMAPS,"./bed2rMAPS/RI_Events_loc_for_rMAPS.txt",sep="\t",quote=F)

