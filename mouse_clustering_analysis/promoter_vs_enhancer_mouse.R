rm(list=ls())
source('mainFunctions_sub.R')
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=28,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=28,face="bold"),
                                 axis.text.x=element_text(size=24),
                                 axis.text.y=element_text(size=24))
# Promoter vs enhancer by tissue ------------------------------------------
#max dMML vs dNME in promoter vs enhancer
enhancer_bin=readRDS(bin_enhancer_rds)
UC_in_dMML_dNME=readRDS(UC_merge_max_loc_file)
UC_in_dMML_dNME=lapply(UC_in_dMML_dNME,function(x) as.data.table(x,keep.rownames=T))
csv_out=data.table()
dMML_in_out=data.table()
dNME_in_out=data.table()
UC_in_out=data.table()
for(fn in dir(dir_out_cluster01)){
  tissue_in=gsub('.csv','',fn)
  if(tissue_in !="NT"){
    csv_in=fread(paste0(dir_out_cluster01,fn))
    olap=findOverlaps(convert_GR(csv_in$regions),enhancer_bin)
    csv_in[,enhancer:=FALSE]
    csv_in[queryHits(olap),enhancer:=TRUE]
    csv_in$enhancer_gene="NA"
    csv_in[queryHits(olap),enhancer_gene:=enhancer_bin$`Target Gene`[subjectHits(olap)]]
    csv_in[,tissue:=tissue_in]
    csv_in[,states:="NA"]
    csv_in[enhancer==TRUE,states:="enhancers"]
    csv_in[abs(distance)<=2000,states:="promoters"]
    csv_out=rbind(csv_out,csv_in)
     dMML_in_out=rbind(dMML_in_out,diff_stat_extraction(UC_in_dMML_dNME[[tissue_in]],'dMML',csv_in))
    dNME_in_out=rbind(dNME_in_out,diff_stat_extraction(UC_in_dMML_dNME[[tissue_in]],'dNME',csv_in))
    UC_in_out=rbind(UC_in_out,diff_stat_extraction(UC_in_dMML_dNME[[tissue_in]],'UC',csv_in))
  }
}
#Boxplot of dMML and dNME

pdf(paste0(figure_path,'promoter_enhancer_dMML_dNME.pdf'))
print(ggplot(dMML_in_out[states%in%c("enhancers","promoters")],aes(x=states,y=dMML ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dMML")+ylim(c(0,0.5))+theme_glob)
print(ggplot(dNME_in_out[states%in%c("enhancers","promoters")],aes(x=states,y=dNME ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dNME")+ylim(c(0,0.5))+theme_glob)
print(ggplot(UC_in_out[states%in%c("enhancers","promoters")],aes(x=states,y=UC ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("UC")+ylim(c(0,0.5))+theme_glob)
dev.off()

#boxplot compare max dNME and dMML
#maximum value
pdf(paste0(figure_path,'promoter_enhancer_dMML_dNME_tissue_max_pair.pdf'))
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dMML_max_pair,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dMML")+theme_glob)
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dNME_max_pair,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dNME")+theme_glob)
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=UC_max_pair,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("UC")+theme_glob)
dev.off()
#at max UC
pdf(paste0(figure_path,'promoter_enhancer_dMML_dNME_tissue_dMML_max_UC.pdf'))
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dMML_max_UC_pair,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dMML")+theme_glob)
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dNME_max_UC_pair ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dNME")+theme_glob)
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=UC_max_pair ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("UC")+theme_glob)
dev.off()

pdf(paste0(figure_path,'promoter_enhancer_dMML_dNME_states_max_pair.pdf'))
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=states,y=dMML_max_pair ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dMML")+theme_glob)
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=states,y=dNME_max_pair ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("dNME")+theme_glob)
print(ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=states,y=UC_max_pair ,fill=states))+geom_boxplot(outlier.shape = NA)+ylab("UC")+theme_glob)
dev.off()
