rm(list=ls())
gc()
library(stringr)
library(data.table)
event_infor = read.table("../../MISO_PSI/merge_PSI/events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",header=T,row.names=1,check.names=F)
isoform_infor <- read.table("../../ref_genome/gencode.vM22.transcript.infor.rm.version.tsv",sep="\t",header=T,row.names=1)
event_in_protein <- apply(event_infor, 1, function(x){
  paste0(isoform_infor[unlist(str_split(x['transcript_id(in)'],';')),"protein_id"],collapse = ";")})
event_ex_protein <- apply(event_infor, 1, function(x){
  paste0(isoform_infor[unlist(str_split(x['transcript_id(ex)'],';')),"protein_id"],collapse = ";")})
event_infor$`protein_id(in)` <- event_in_protein
event_infor$`protein_id(ex)` <- event_ex_protein

EC_T1_included = read.table("../../../NewAnalysis/MISO_PSI/Anchor/EC_T1/EC_T1_signature_included_events_with_UTR_Srsf2_site.txt",header=T,row.names=1,sep="\t")
EC_T1_included_events = rep("EC_T1_included",nrow(EC_T1_included))
names(EC_T1_included_events) = row.names(EC_T1_included)

fish_candi = read.delim("../../../NewAnalysis/MISO_PSI/FISH_confirm/FISH_supplemental.txt",sep="\t",header=T)
new_fish_candi = read.table("../../../NewAnalysis/MISO_PSI/FISH_confirm/T1_included_new_8_candidates.txt",sep="\t",header=T)
all_fish_candi = union(as.character(fish_candi$event_name),
as.character(new_fish_candi$event_name))
fish_candi_tab = rep("FISH",length(all_fish_candi))
names(fish_candi_tab) = all_fish_candi

events_modality = read.table("./EventModality/7_stage_event_modality.txt",sep="\t",header=T,row.names=1)
events_modality["chr5:100375208:100375331:-@chr5:100373425:100373715:-@chr5:100368294:100368352:-",]
EC_T1_modality = events_modality[,c("AEC","HEC","T1_pre_HSC")]
#EC_T1_modality = EC_T1_modality[apply(EC_T1_modality,1,function(x){length(x[x!=""])<3}),]

dir.create("AEC_HEC_T1")
#------AEC HEC included
AEC_HEC_in = EC_T1_modality[EC_T1_modality$HEC=="included",]
AEC_HEC_in = AEC_HEC_in[AEC_HEC_in$AEC!="included",]
AEC_HEC_in_out = data.frame(event_name = row.names(AEC_HEC_in), events_modality[row.names(AEC_HEC_in),], 
event_infor[row.names(AEC_HEC_in),],check.names=F)

AEC_HEC_in_out$is_T1_sig = EC_T1_included_events[row.names(AEC_HEC_in_out)]
AEC_HEC_in_out$FISH = fish_candi_tab[row.names(AEC_HEC_in_out)]
fwrite(as.data.table(AEC_HEC_in_out),"./AEC_HEC_T1/AEC_HEC_included_events.txt",sep="\t",quote=F,row.names=F)
AEC_HEC_T1_in_out = AEC_HEC_in_out[AEC_HEC_in_out$T1_pre_HSC=="included",]
fwrite(as.data.table(AEC_HEC_T1_in_out),"./AEC_HEC_T1/AEC_other_HEC_T1_included_events.txt",sep="\t",quote=F,row.names=F)

#------AEC T1 included
AEC_T1_in = EC_T1_modality[EC_T1_modality$T1_pre_HSC=="included",]
AEC_T1_in = AEC_T1_in[AEC_T1_in$AEC!="included",]
AEC_T1_in_out = data.frame(AEC_T1_in, event_infor[row.names(AEC_T1_in),],check.names=F)
AEC_T1_in_out$is_T1_sig = EC_T1_included_events[row.names(AEC_T1_in_out)]
AEC_T1_in_out$FISH = fish_candi_tab[row.names(AEC_T1_in_out)]
fwrite(as.data.table(AEC_T1_in_out),"./AEC_HEC_T1/AEC_T1_included_events.txt",sep="\t",quote=F,row.names=F)

#------HEC T1 included
HEC_T1_in = EC_T1_modality[EC_T1_modality$T1_pre_HSC=="included",]
HEC_T1_in = HEC_T1_in[HEC_T1_in$HEC!="included",]
HEC_T1_in_out = data.frame(event_name = row.names(HEC_T1_in), events_modality[row.names(HEC_T1_in),],
event_infor[row.names(HEC_T1_in),],check.names=F)
HEC_T1_in_out$is_T1_sig = EC_T1_included_events[row.names(HEC_T1_in_out)]
HEC_T1_in_out$FISH = fish_candi_tab[row.names(HEC_T1_in_out)]
fwrite(as.data.table(HEC_T1_in_out),"./AEC_HEC_T1/HEC_T1_included_events.txt",sep="\t",quote=F,row.names=F)
AEC_HEC_other_T1_in = HEC_T1_in_out[HEC_T1_in_out$AEC!="included",]
fwrite(as.data.table(AEC_HEC_other_T1_in),"./AEC_HEC_T1/T1_included_AEC_HEC_not_included_events.txt",sep="\t",quote=F,row.names=F)

#------HEC included and T1 included
HEC_in = EC_T1_modality[EC_T1_modality$HEC=="included",]
HEC_in_T1_not = HEC_in[HEC_in$T1_pre_HSC!="included",]
HEC_in_T1_not_out = data.frame(HEC_in_T1_not, event_infor[row.names(HEC_in_T1_not),],check.names=F)
HEC_in_T1_not_out$is_T1_sig = EC_T1_included_events[row.names(HEC_in_T1_not)]
HEC_in_T1_not_out$FISH = fish_candi_tab[row.names(HEC_in_T1_not)]
fwrite(as.data.table(HEC_in_T1_not_out),"./AEC_HEC_T1/HEC_included_T1_not_included_events.txt",sep="\t",quote=F)









