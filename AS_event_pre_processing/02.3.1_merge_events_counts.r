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
psi_tab = all_event[,c("event_name","assigned_counts")]
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
outPath = "./merged_counts"

all_stage = c("AEC","HEC","T1_pre_HSC","T2_pre_HSC","E12","E14","Adult_HSC")
#running merge PSI steps, it takes few minutes.
all_PSI = getStagePSI(stage = all_stage[1])
for(i in all_stage[-1]){
psi = getStagePSI(stage = i)
all_PSI = merge(x = all_PSI, y = psi, by = "event_name", all=TRUE)
}

fwrite(all_PSI, file.path(outPath,'all_event_counts_Tab.txt'),sep="\t",quote=F,row.names=F)

all_PSI_frm = fread(file.path(outPath,'all_PSI_Tab.txt'),sep="\t",header=T,data.table=F)
all_PSI = all_PSI_frm[,-1]
all_PSI_num = apply(all_PSI,2,function(x){length(x[which(!is.na(x))])})
write.table(all_PSI_num, 'all_sample_event_number.txt',sep="\t",quote=F)


all_counts_frm <- all_PSI
all_counts_tab <- all_counts_frm[,-1]
row.names(all_counts_tab) <- all_counts_frm$event_name

sum_mat <- matrix(NA,nrow = nrow(all_counts_tab),ncol =  ncol(all_counts_tab))
for(i in 1:nrow(all_counts_tab)){
for(j in 1:ncol(all_counts_tab)){
if(!is.na(all_counts_tab[i,j])){
y <- unlist(strsplit(all_counts_tab[i,j],','))
z <- sapply(y,function(x){unlist(strsplit(x,':'))[2]})
sum_mat[i,j] <- sum(as.numeric(z))
}
}
}

dimnames(sum_mat) <- dimnames(all_counts_tab)

sum_mat_median <- apply(sum_mat,1,function(x){mean(x,na.rm=T)})
plot_frm <- data.frame(counts = sum_mat_median)

sum_mat_all <- as.vector(sum_mat)
sum_mat_all <- na.omit(sum_mat_all)
plot_frm <- data.frame(counts = sum_mat_all)


library(ggplot2)
line_size = 0.5
text_size = 8
mytheme <- theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, hjust = 0.5, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')


p <- ggplot(data = plot_frm, aes(log2(counts+1)))+
geom_density()+
#geom_histogram(bins=60,fill="#f48fb1")+
geom_vline(aes(xintercept = log2(6)),color="red", linetype="dotted")+
xlab("uniquely mapped reads (X 10^6)")+
#ylab("number of AS evvents")+
#scale_x_continuous(limits=c(0,1000))+
#scale_y_continuous(limits=c(0,9000),breaks = seq(0,8000,by=2000))+
mytheme


p





