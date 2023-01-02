rm(list=ls())
library(xlsx)
GO_in=as.data.table(read.xlsx('../downstream/output/graphs_tables/regulatory_non_regulatory_GOrilla.xlsx',1))

for(i in 1:nrow(GO_in)){
  if(!is.na(GO_in$Description[i])){GO_term=GO_in$Description[i]}
  else{GO_in$Description[i]=GO_term}
  
  
}
GO_in=GO_in[,list(P.value=max(P.value,na.rm=T),
            FDR.q.value=max(FDR.q.value,na.rm=T),
            Enrichment..N..B..n..b.=max(as.numeric(gsub('\\(.*','',Enrichment..N..B..n..b.)),na.rm=T),
            Genes=gsub('\\[-] Hide genes,','',paste(unique(gsub(' - .*','',Genes)),collapse =','))),by=Description]
top50_gene=read.xlsx2('../downstream/output/human_analysis/NME_motif/NME_regulaotry_Ken.xlsx',1,startRow = 2)

GO_in$Genes=unlist(lapply(strsplit(GO_in$Genes,','),function(x) paste(x[x%in%top50_gene[1:50,]$Transcription.facotrs],collapse = ',')))

write.csv(GO_in,'../downstream/output/graphs_tables/regulatory_non_regulatory_GOrilla_processed.csv')
