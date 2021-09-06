rm(list=ls())
gc()
library(data.table)
assignModality <- function(stage, log2bf_cutoff, inpath = anchor_out_path, outpath = outPath){
modality_tab = read.csv(paste(inpath,paste(stage,'_PSI_Filter_5_Modality.csv',sep=''),sep='/'),header=T,row.names=1)
modality_tab[is.na(modality_tab)]=0
modality_frm = data.frame(event_name = row.names(modality_tab),modality_tab,check.rows=T)
modality_melt = melt(modality_frm, id.vars = c("event_name"), variable.name="modality",value.name="log2BF")
modality_melt_frm = as.data.frame(modality_melt)
modality_mat = modality_melt_frm[modality_melt_frm$log2BF >= log2bf_cutoff,]

#------get bimodal events.
bimodal_events = modality_mat[modality_mat$modality=="bimodal",]
row.names(bimodal_events) = bimodal_events$event_name
#------get excluded events. 
excluded_events = modality_mat[modality_mat$modality=="excluded",]
row.names(excluded_events) = excluded_events$event_name
#------get included events.
included_events = modality_mat[modality_mat$modality=="included",]
row.names(included_events) = included_events$event_name
#------get middle events.
middle_events = modality_mat[modality_mat$modality=="middle",]
row.names(middle_events) = middle_events$event_name
#------get uncategorized events.
uncategorized_events = modality_mat[modality_mat$modality=="uncategorized",]
row.names(uncategorized_events)  = uncategorized_events$event_name
#------merge 4 modality. 
four_modality_events = unique(c(row.names(bimodal_events), row.names(excluded_events),
row.names(included_events), row.names(middle_events)))
#------get multimodal events.
multimodal_event_name = setdiff(row.names(uncategorized_events), four_modality_events)
multimodal_event_frm = data.frame(event_name = multimodal_event_name,modality = rep('multimodal',length(multimodal_event_name)))
four_modality_frm = modality_tab[four_modality_events,1:(ncol(modality_tab)-1)]
#------get five types of modality.
events_four_modality = apply(four_modality_frm, 1, function(x){names(which(x == max(x)))})
events_four_modality_frm = data.frame(event_name = names(events_four_modality),modality = events_four_modality)
all_event_frm = rbind(events_four_modality_frm, multimodal_event_frm)
modality_dir = paste(outpath,'Events2Modality',sep='/')
dir.create(modality_dir, recursive = T, showWarnings=F)
write.table(all_event_frm,paste(modality_dir,paste(stage,"_all_events_modality.txt",sep=""),sep='/'),sep="\t",quote=F,row.names=F)
#------write all events modality Bayes Factor table.
four_modality_BF_tab = t(apply(four_modality_frm,1,function(x){c(paste(names(x[which(x >= log2bf_cutoff)]),collapse=';'),
paste(x[which(x >= log2bf_cutoff)],collapse=';'))}))
four_modality_BF_frm = data.frame(event_name = row.names(four_modality_BF_tab), four_modality_BF_tab, stringsAsFactors = F)
colnames(four_modality_BF_frm) = c("event_name","modality","log2BF")
multimodal_modality_BF_frm = data.frame(event_name = multimodal_event_name, modality = rep("multimodal", length(multimodal_event_name)),
log2BF = modality_tab[multimodal_event_name,"uncategorized"], stringsAsFactors = F)
all_event_BF_frm = rbind(four_modality_BF_frm, multimodal_modality_BF_frm)
BF_dir = paste(outpath,'EventsModalityBF',sep='/')
dir.create(BF_dir, recursive = T, showWarnings=F)
write.table(all_event_BF_frm,paste(BF_dir,paste(stage,"_all_events_modality_bayes_factor.txt",sep=""),sep='/'),sep="\t",quote=F,row.names=F)

#return(all_event_frm)
}

anchor_out_path = "./AnchorOut"
outPath = "./EventModality"
dir.create(outPath)
stages = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
for(i in stages){
assignModality(stage = i,log2bf_cutoff=1)
}




