library(data.table)
event_infor = fread("./events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",
header=T,data.table=F,check.names=F,stringsAsFactors=F)
row.names(event_infor) = event_infor$event_name

PSI_tab = fread("./merged_PSI/all_PSI_Tab.txt",sep="\t",header=T,data.table=F,check.names=F,stringsAsFactors=F)
row.names(PSI_tab) = PSI_tab$event_name

other_events = setdiff(row.names(PSI_tab), row.names(event_infor))
left_events = PSI_tab[other_events,]
left_events_num = apply(left_events,1,function(x){length(x[!is.na(x)])})

