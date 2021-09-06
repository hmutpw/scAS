library(data.table)
dir = "./events2gene/events2transcript"
list.files(dir)

SE_tab = fread(file.path(dir,"SE","SE_events_to_transcript_id.txt"),sep="\t",header=T,data.table=F,stringsAsFactors=F)
SE_tab$event_type = rep("SE",nrow(SE_tab))
row.names(SE_tab) = SE_tab$event_name
A3SS_tab = fread(file.path(dir,"A3SS","A3SS_events_to_transcript_id.txt"),sep="\t",header=T,data.table=F,stringsAsFactors=F)
A3SS_tab$event_type = rep("A3SS",nrow(A3SS_tab))
row.names(A3SS_tab) = A3SS_tab$event_name
A5SS_tab = fread(file.path(dir,"A5SS","A5SS_events_to_transcript_id.txt"),sep="\t",header=T,data.table=F,stringsAsFactors=F)
A5SS_tab$event_type = rep("A5SS",nrow(A5SS_tab))
row.names(A5SS_tab) = A5SS_tab$event_name
MXE_tab = fread(file.path(dir,"MXE","MXE_events_to_transcript_id.txt"),sep="\t",header=T,data.table=F,stringsAsFactors=F)
MXE_tab$event_type = rep("MXE",nrow(MXE_tab))
row.names(MXE_tab) = MXE_tab$event_name
RI_tab = fread(file.path(dir,"RI","RI_events_to_transcript_id.txt"),sep="\t",header=T,data.table=F,stringsAsFactors=F)
RI_tab$event_type = rep("RI",nrow(RI_tab))
row.names(RI_tab) = RI_tab$event_name

all_event_tab = rbind(SE_tab,A3SS_tab,A5SS_tab,MXE_tab,RI_tab)


#get isoform infor.
isoInfor = fread("../../ref_genome/gencode.vM22.transcript.infor.tsv",sep="\t",header=T,data.table=F,stringsAsFactors=F)
row.names(isoInfor) = isoInfor$transcript_id
geneInfor = fread("../../ref_genome/gencode.vM22.gene.infor.tsv",sep="\t",header=T,data.table=F,stringsAsFactors=F)
row.names(geneInfor) = geneInfor$gene_id

#these 3 steps will take few minutes
all_event_tab$in_transcript_name = apply(all_event_tab,1,function(x){
paste0(isoInfor[unlist(strsplit(x['included_transcript'],split=';',fixed=T)),"transcript_name"],collapse=";")})
all_event_tab$ex_transcript_name = apply(all_event_tab,1,function(x){
paste0(isoInfor[unlist(strsplit(x['excluded_transcript'],split=';',fixed=T)),"transcript_name"],collapse=";")})
all_event_tab$gene_name = apply(all_event_tab,1,function(x){
geneInfor[unlist(strsplit(x['gene_id'],split=';',fixed=T)),"gene_name"]})
#remove version number.
all_event_tab$gene_id = substr(all_event_tab$gene_id,1,18)
all_event_tab$included_transcript = apply(all_event_tab,1,function(x){
paste0(substr(unlist(strsplit(x['included_transcript'],split=';',fixed=T)),1,18),collapse=";")})
all_event_tab$excluded_transcript = apply(all_event_tab,1,function(x){
paste0(substr(unlist(strsplit(x['excluded_transcript'],split=';',fixed=T)),1,18),collapse=";")})

all_event_tab_out = all_event_tab[,c("event_name","event_type","exon_num",
"gene_id","gene_name","included_transcript","in_transcript_name",
"excluded_transcript","ex_transcript_name")]
colnames(all_event_tab_out) = c("event_name","event_type","exon_num",
"gene_id","gene_name",'transcript_id(in)','transcript_name(in)',
'transcript_id(ex)','transcript_name(ex)')
fwrite(all_event_tab_out,"./events2gene/gencode.vM22.MISO.Events2gene.txt",sep="\t",quote=F,row.names=F)





