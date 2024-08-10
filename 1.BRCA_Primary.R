library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(infercnv);library(CaSpER);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
my36colors <-c('#3638a9','#e8e967','#f39562','#c7eedb','#c496d4','#7e7913','#00223d','#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1",'red','red','red')

mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1)),
        axis.text.x = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))

set.seed(1234)

main.path='/data4/huanggy/Fusion/GSE180286_BRCA_Onco2021/'
save.data='/data4/huanggy/Fusion/GSE180286_BRCA_Onco2021/data/'
save.pic='/data4/huanggy/Fusion/GSE180286_BRCA_Onco2021/pic/'

#---marker-----
major.marker.list=list(Tcell=c('CD2','CD3D','CD3E','CD3G','CD247','CD7','CD8A','CD8B','CD4'),
                       Bcell=c('CD19','MS4A1','CD79A','CD79B','BANK1'),Plasma=c('MZB1','TNFRSF17','SDC1','XBP1'),
                       Myeloid=c('CD14','FCGR3A','C1QA','C1QB','MARCO','CD163','FCGR2A','CD1E','CD1A','CD68','S100A8','S100A9','S100A12','LILRA4','IL3RA','VCAN','CLEC9A','XCR1'),
                       Mast=c('MS4A2','CPA3','TPSAB1','TPSB2'),Fibroblast=c('COL1A2','COL3A1','COL1A1','COL6A2','DCN','COL6A3','COL6A1','COL5A1','MMP2','ACTA2','PDGFRA','CNN1'),
                       Endothelial=c('ENG','CLDN5','PECAM1','CDH5','VWF','FLT1'),Epithelial=c('EPCAM','KRT19','KRT18','KRT8','KRT7')
)
major.marker.list=data.frame(unlist(major.marker.list) ,celltype= names(unlist(major.marker.list)) %>% str_remove('\\d+') , row.names = as.vector(unlist(major.marker.list)) )



#-----Comnine data------
Patients.use=list.files(paste0(main.path,'/raw/'))

raw.counts=list()
for(i in Patients.use){
  tmp=Read10X(paste0(main.path,'/raw/',i,'/','/filtered_feature_bc_matrix'))
  tmp@Dimnames[[2]]=paste0('BRCA_',i,'_Ds0Onco20210_',tmp@Dimnames[[2]]) %>% str_extract('\\w+')
  tmp@Dimnames[[2]]<- str_replace_all(tmp@Dimnames[[2]],'-','_')
  raw.counts[[i]]=tmp
}

#-----delete mito------
counts=mclapply(1:length(raw.counts),function(i){
  tmp <- CreateSeuratObject(counts =raw.counts[[i]])
  tmp[["percent.mt"]] = PercentageFeatureSet(tmp, pattern = "^[Mm]T-")
  tmp[["RP.mt"]] = PercentageFeatureSet(tmp, pattern = "^RP[LS]\\w")
  tmp <- subset(tmp, subset = (nFeature_RNA > 300 & nFeature_RNA< 5000 & percent.mt < 20 & RP.mt < 50))
  tmp=tmp@assays$RNA@counts
  return(tmp)
},mc.cores = getOption('mc.cores',10))

all.count.movedoublet=list()
for(i in 1:length(counts)){
  tmp <- CreateSeuratObject(counts[[i]]) %>% NormalizeData(verbose =F) %>% ScaleData(verbose =F) %>%  FindVariableFeatures(selection.method = "vst", nfeatures = 2000,verbose =F) %>% RunPCA(verbose =F) %>% FindNeighbors(dims =1:10,verbose =F) %>% FindClusters(resolution =0.5,verbose =F)%>% RunUMAP( dims = 1:10,verbose =F)
  homotypic.prop <- modelHomotypic(tmp@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.05*(ncol(tmp)) ) ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:10, pN = 0.2, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  tmp=tmp@meta.data
  colnames(tmp)[6]='PANN'
  colnames(tmp)[7]='DF'
  all.count.movedoublet[[i]]=tmp
}

for(i in 1:length(counts) ){
  choose.barcode=all.count.movedoublet[[i]][all.count.movedoublet[[i]]$DF=='Singlet',] %>% rownames()
  counts[[i]]=counts[[i]][,choose.barcode]
}

saveRDS(counts,paste0(save.data,'Alldata.removedoublet.rds'))


#------cluster step1--------

counts=Reduce(cbind,counts)
counts <- CreateSeuratObject(counts =counts)
counts[["percent.mt"]] = PercentageFeatureSet(counts, pattern = "^[Mm]T-")
counts[["percent.RPL"]] = PercentageFeatureSet(counts, pattern = "^RP[SL].+")

counts <- NormalizeData(object = counts, normalization.method = "LogNormalize", scale.factor = 10000)
counts <- FindVariableFeatures(object = counts, selection.method = "vst", nfeatures = 2000)
counts <- ScaleData(counts, vars.to.regress = c("nCount_RNA"))
counts <- RunPCA(object =counts, features = VariableFeatures(object = counts))


message=rownames(counts@meta.data) %>% str_split_fixed('_',n=5) 
counts@meta.data[,'Cancer']=message[,1];counts@meta.data[,'Patient']=message[,2]
counts@meta.data[,'Tissue']=message[,3];counts@meta.data[,'Dataset']=message[,4]


counts <- RunHarmony(counts, group.by.vars = c("Patient"), dims.use = 1:33, verbose = T)
counts <- FindNeighbors(object = counts, dims = 1:35, reduction = "harmony")
counts <- FindClusters(object = counts, resolution = 0.5)
counts<- RunUMAP(object = counts, dims = 1:30, reduction = "harmony")

saveRDS(counts,paste0(save.data,'Alldata.firstcluster.rds'))



#----------define subtype-------------


markers=wilcoxauc(counts) %>% { .[.$logFC>0.8 & .$pval<0.05,]  }
markers=markers[markers$feature %in% rownames(major.marker.list), ] %>% {  dplyr::mutate( . , Sug.type=major.marker.list[.$feature,'celltype'] )  } 

define.cluster=markers %>% data.table() %$% .[,.N,.(group,Sug.type)] %>% data.frame()
final.cluster=list()
for(i in unique(define.cluster$group)  ) {
  now.define.cluster=define.cluster[define.cluster$group==i,]
  
  if( (nrow(now.define.cluster)==1 & now.define.cluster[,'N']>1)  ){final.cluster[[i]]=c(i,now.define.cluster$Sug.type) 
  }else{
    if( max(now.define.cluster$N) - sum(setdiff(now.define.cluster$N, max(now.define.cluster$N)  )) >1  ) {
      final.cluster[[i]]= c(i, now.define.cluster[now.define.cluster$N==max(now.define.cluster$N),'Sug.type'] )
    }
  }
  
  # if(length(final.cluster[[i]])>2){
  #   final.cluster=final.cluster[names(final.cluster)!='15']
  # }
  
}
final.cluster=data.frame(final.cluster) %>% t() %>% data.frame() %>% { rownames(.) = .[,1]; colnames(.)= c('seurat_clusters','Sug.type') ; .}
final.cluster=rbind(final.cluster,data.frame(seurat_clusters= setdiff(counts$seurat_clusters , final.cluster[,1] ) ,Sug.type =NA) %>% { rownames(.)=.[,1];.} )

counts$recluster=NA
counts$recluster=final.cluster[as.character(counts$seurat_clusters),'Sug.type']

#----report1-----
Print.report(counts,c('CD3D','CD19','TNFRSF17','CD14','CPA3','PECAM1','COL6A1','EPCAM'), 'All.cell.Report' )
print(paste0('现阶段尚未被注释的群:',paste(final.cluster[is.na(final.cluster$Sug.type),1],collapse='-')))


final.cluster=data.frame(seurat_clusters=final.cluster[is.na(final.cluster$Sug.type),'seurat_clusters'], Sug.type=c('Epithelial','Tcell','Myeloid','Myeloid','Epithelial','Tcell','Epithelial','Myeloid','Epithelial')) %>% { rbind(.,final.cluster)} %>%   {.[!is.na(.[,2]),] } %>% { rownames(.)=.$seurat_clusters ; .}
counts$recluster=final.cluster[as.character(counts$seurat_clusters),'Sug.type']

saveRDS(counts,paste0(save.data,'Alldata.twocluster.rds'))



#-----cluster Epithelial----
count.onetype=subset(counts,subset=recluster=='Epithelial')
count.onetype=Delete.recluster(count.onetype,FindNei.dim=35)
Print.report(count.onetype,c('CD3D','CD3E','CD19','CD14','CD1C','LILRA4','IL3RA','CLEC9A','PECAM1','EPCAM','KRT18','KRT19','KRT8','KRT7','MUC6','TFF3','COL6A1','VIM','CDH2','TP63'), 'Epi.report' ,show.legend=T)

#markers=wilcoxauc(count.onetype,group_by = 'seurat_clusters') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}

count.onetype=subset(count.onetype,subset=seurat_clusters %in% c(0:3,5:8,10:11,14,15) )
count.onetype=Delete.recluster(count.onetype,FindNei.dim=35)
Print.report(count.onetype,c('CD3D','CD3E','CD19','CD14','CD1C','LILRA4','IL3RA','CLEC9A','PECAM1','EPCAM','KRT18','KRT19','KRT8','KRT7','MUC6','TFF3','COL6A1','VIM','CDH2','TP63'), 'Epi.report2' ,show.legend=T)

count.onetype=subset(count.onetype,subset=seurat_clusters %in% c(0:8,10:13) )
count.onetype=Delete.recluster(count.onetype,FindNei.dim=35)
Print.report(count.onetype,c('CD3D','CD3E','CD19','CD14','CD1C','LILRA4','IL3RA','CLEC9A','PECAM1','EPCAM','KRT18','KRT19','KRT8','KRT7','MUC6','TFF3','COL6A1','VIM','CDH2','TP63'), 'Epi.report3' ,show.legend=T)

count.onetype=subset(count.onetype,subset=seurat_clusters %in% c(0:8,10,12,13) )
count.onetype=Delete.recluster(count.onetype,FindNei.dim=35)
Print.report(count.onetype,c('CD3D','CD3E','CD19','CD14','CD1C','LILRA4','IL3RA','CLEC9A','PECAM1','EPCAM','KRT18','KRT19','KRT8','KRT7','MUC6','TFF3','COL6A1','VIM','CDH2','TP63'), 'Epi.report4' ,show.legend=T)

count.onetype$recluster=paste0('Epi_C',count.onetype$seurat_clusters)

gene.signature=data.frame( unique(count.onetype$recluster %>% str_extract('C\\d+') ), gene.choose=c('ERBB2','STMN1','STMN1','KRT7','LPR1B','KCNJ3','MMP7','JUN','PTPN1','CP','ANKRD30A')) %>% { row.names(.)=.[,1];.}
count.onetype$recluster =paste0(count.onetype$recluster %>% str_split_fixed('_',n=2) %>% { .[,1]}, '_',gene.signature[count.onetype$recluster %>% str_extract('C\\d+'),'gene.choose'] )
Print.report(count.onetype,c('EPCAM','KRT18','SFTPB','AGER','SCGB1A1','FOXJ1','KRT14','ASCL1','VIM'), 'Epi.report9',show.legend=T )

markers=wilcoxauc(count.onetype,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
top2heatmap=top2heatmap[!duplicated(top2heatmap$feature),]
marker_heatmap2(top2heatmap,count.onetype,save.pic,'recluster','heatmap.Epi4', levels(as.factor(top2heatmap$group)) ,height=20 )
dev.off()

saveRDS(count.onetype,paste0(save.data,'Epi.rds'))

#----infercnv-----
system(paste0('/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ',main.path,'program/CNV.method1.All.R ',main.path,'/ '),wait=F,show.output.on.console=F)


#----copykat------
system(paste0('/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ',main.path,'program/CNV.copy.All.R ',main.path,'/ '),wait=F,show.output.on.console=F)
patient= count.onetype$Patient %>% table() %>% data.frame() %>% {.[.$Freq>30,'.']} %>% as.character()
for(i in unique(count.onetype@meta.data[count.onetype$Tissue %in% c('Tumor','LN'),'Patient'])[31:41] ){
  system(paste0('/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ',main.path,'program/CNV.method1.patient.R ',main.path,'/ ',i),wait=F,show.output.on.console=F)
}

for(i in unique(count.onetype@meta.data[count.onetype$Tissue %in% c('Tumor','LN'),'Patient'])[31:41] ){
  system(paste0('/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ',main.path,'program/CNV.copy.patient.R ',main.path,'/ ',i),wait=F,show.output.on.console=F)
}


copykat.all=read.table(paste0(save.data,'/cnv/copy/All/All_copykat_prediction.txt'),header=T,row.names = 1)
count.onetype@meta.data[rownames(copykat.all),'copykat.all']=copykat.all$copykat.pred
count.onetype@meta.data[count.onetype$Tissue %in% c('Normal','Blood','Ptumor'),'copykat.all']='diploid'

for(i in list.files(paste0(save.data,'/cnv/copy/')) %>%{setdiff(.,'All')} ){
  copykat.one=read.table(paste0(save.data,'/cnv/copy/',i,'/',i,'_copykat_prediction.txt'),header=T,row.names = 1)
  count.onetype@meta.data[rownames(copykat.one),'copykat.patient']=copykat.one$copykat.pred
}
count.onetype@meta.data[count.onetype$Tissue %in% c('Normal','Blood','Ptumor'),'copykat.patient']='diploid'


infercnv.all=read.table(paste0(save.data,'/cnv/infercnv/All/map_metadata_from_infercnv.txt'),header=T,row.names = 1)
count.onetype@meta.data[rownames(infercnv.all),paste0('has_cnv_chr',c(seq(1,22),'X','Y'))]=infercnv.all[,paste0('has_cnv_chr',c(seq(1,22),'X','Y'))]
count.onetype@meta.data[rownames(infercnv.all),'infercnv.all']=apply(infercnv.all[,paste0('has_cnv_chr',c(seq(1,22),'X','Y'))],1,sum)


for(i in list.files(paste0(save.data,'/cnv/infercnv/')) %>%{setdiff(.,'All')} ){
  infercnv.one=read.table(paste0(save.data,'/cnv/infercnv/',i,'/','map_metadata_from_infercnv.txt'),header=T,row.names = 1)
  count.onetype@meta.data[rownames(infercnv.one),'infercnv.patient']=apply(infercnv.one[,paste0('has_cnv_chr',c(seq(1,22),'X','Y'))],1,sum)
}

saveRDS(count.onetype,paste0(save.data,'Epi.rds'))

