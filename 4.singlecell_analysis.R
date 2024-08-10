.libPaths('/data3/R_package2')
dyn.load('/data3/anaconda3/lib/libhdf5_hl.so.100')
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(copykat);#library(CaSpER);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library(presto);library(harmony);library(ggrepel)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
library(homologene)
library(KEGGREST)
library(pathview) 
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(infercnv);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(Startrac);library(scRepertoire)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
library(ROGUE)
library(SingleCellExperiment)
library(smoother)
library(plotrix)
library(ggpubr)
library(RColorBrewer)
library(Seurat)
library(velocyto.R)
library(magrittr)
library(stringr)
library(SeuratDisk)
library(ggalluvial)
library(ggforce)
library(SeuratWrappers)
library(org.Mm.eg.db)
library(AUCell)
library(ape)
library(loomR)
library(Hmisc)
library(irGSEA)
library(CytoTRACE)
library(genomicInstability)
library(Chord)
library(scMetabolism)
library(projectLSI)
library(ggridges)
#BiocManager::install("maftools")
library(maftools)
library(CellChat)
library(ConsensusClusterPlus)
library(loomR)
library(SCopeLoomR)
library(SCENIC)
library(plot1cell)
library(URD)
library(monocle3)
library(tidyr)
library(clusterProfiler)
library(GseaVis)
library(org.Mm.eg.db)

library(enrichplot)
library(tidyverse)
library(ggstatsplot)

my36colors <-c('#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1",'red','blue','green','red','red')
my36colors1<-c('#594C57',"#8F6D5F", '#EAE5E3','#DB9C5E',"#420B2F","#C25160","#CC5D20",'#2E59A7','#80A492','#32788A','#E1D279','#A64036')

my36colors2 <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"))  

my36colors3 <- c('#EAE5E3','#594C57',"#8F6D5F","#420B2F","#C25160","#CC5D20",'#2E59A7','#80A492','#32788A','#E1D279','#A64036','#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
                 "mediumvioletred","gold", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1",'red','blue','green',brewer.pal(9, "Set1"),'grey')  

set.seed(1234)

main.path='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/'
save.data='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/'
save.pic='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/'



#----infercnv-----
.libPaths('/data1/liuxf/R_package4')
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder)
library(presto);library(harmony);library(ggrepel);library(copykat)
library(infercnv)
#Args <- commandArgs(T)

#Args <- c("/data10/liuxf/PRAD/All1/Epi/",'D1_P10')
#save.data=Args[1]
save.data='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/'
system(paste0('mkdir -p ',save.data,'cnv/infercnv/TREM2'))

counts <- readRDS("/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Alldata.twocluster.rds")
Fib <- readRDS("/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Fib.rds")
Epi <- readRDS("/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Epi.rds"


Epi_count.cnv=counts
#Epi_count.cnv@meta.data[  rownames(Onco2021_Epi) ,'recluster1']=Onco2021_Epi@meta.data[ rownames(Onco2021_Epi) ,'recluster']

set.seed(1000)

choose.barcode= union(rownames(Epi@meta.data[Epi$recluster %in% c('Epi_TREM2'),])  , sample(rownames(Epi_count.cnv@meta.data[Epi_count.cnv$recluster %in% c('Fibroblast'),]),100)) 
Epi_count.cnv=Epi_count.cnv %>% subset(cells=choose.barcode)
Epi_count.cnv@meta.data[rownames(Epi_count.cnv@meta.data),'recluster1'] <- Epi_count.cnv@meta.data$recluster
Epi_count.cnv@meta.data[Epi_count.cnv$recluster %in% c('Fibroblast'),'recluster1']='Normal'


Epi_count.cnv=Epi_count.cnv %>% subset(cells= rownames(Epi_count.cnv@meta.data[!Epi_count.cnv$recluster1 %in% delete.cluster,] ) )
ann=read.csv('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/program/now_gene_pos.txt',header=F,sep='\t')
count=Epi_count.cnv@assays$RNA@counts
intersectgene=intersect(rownames(count),ann$V1)
ann_new=ann[ann$V1 %in% intersectgene,]
count_new=count[ann_new$V1,] %>% data.frame()




meta=data.frame(row.names = colnames(Epi_count.cnv),Epi_count.cnv$recluster1)
ann_new=data.frame(row.names = ann_new$V1,ann_new[2:4])

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_new,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=ann_new,
                                    ref_group_names=c("Normal"),
                                    chr_exclude=c('chrM')) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0(save.data,'cnv/infercnv/TREM2'),  
                             cluster_by_groups=T,
                             num_threads=10,
                             noise_filter=0.2,
                             denoise=T,
                             HMM=T) 

add_to_seurat(seurat_obj = Epi_count.cnv, paste0(save.data,'cnv/infercnv/TREM2'), top_n = 10,bp_tolerance = 2e+06)



#-------genomicInstabilityScore-----

count.one.type.epi=Epi_count.cnv

addDensity <- function(x, col="grey") {
  polygon(c(x$x[1], x$x, x$x[length(x$x)], x$x[1]), c(0, x$y, 0, 0), col=col)
}

genomicInsta=count.one.type.epi
genomicInsta=genomicInsta 

genomicInsta@meta.data[,'barcode']=rownames(genomicInsta@meta.data)
genomicInsta=subset(genomicInsta,cells=  (doBy::sampleBy(formula = ~ recluster,frac = 0.3,data =genomicInsta@meta.data) %>% {.$barcode} )   )
genomicInsta1=genomicInsta
genomicInsta=as.matrix(genomicInsta@assays$RNA@data)
nameschange=bitr(rownames(genomicInsta),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
nameschange=nameschange[!duplicated(nameschange$SYMBOL),]
rownames(nameschange)=nameschange$SYMBOL
genomicInsta=data.frame(nameschange[rownames(genomicInsta),2],genomicInsta)
genomicInsta=genomicInsta[!is.na(genomicInsta$nameschange.rownames.genomicInsta...2.),]
genomicInsta=genomicInsta[!duplicated(genomicInsta$nameschange.rownames.genomicInsta...2.),]
rownames(genomicInsta)=genomicInsta[,1]
genomicInsta=genomicInsta[,-1]

cnv_norm <- inferCNV(as.matrix(genomicInsta), nullmat= as.matrix(genomicInsta)[,colnames(as.matrix(genomicInsta)) %in% rownames(Fib@meta.data)]  )
cnv_norm <- genomicInstabilityScore(cnv_norm, likelihood=TRUE)
giDensityPlot(cnv_norm, ylim=c(0, 2.5))

pos <- order(cnv_norm$gis)
lines(cnv_norm$gis[pos], cnv_norm$gi_likelihood[pos], lwd=2, col="blue")
axis(4, seq(0, 1, length=6), seq(0, 1, .2), col="blue", col.axis="blue")
axis(4, .5, "Relative likelihood", tick=FALSE, line=1.5, col.axis="blue")
pos5 <- which.min((.5-cnv_norm$gi_likelihood)^2)
lines(c(rep(cnv_norm$gis[pos5], 2), max(cnv_norm$gis*1.05)),
      + c(0, rep(cnv_norm$gi_likelihood[pos5], 2)), lty=3, col="blue")



gis_RGS5 <- cnv_norm$gis[genomicInsta1@meta.data[names(cnv_norm$gis),'recluster'] %in% c("Epithelial") ]
gis_RGS51 <- cnv_norm$gis[genomicInsta1@meta.data[names(cnv_norm$gis),'recluster'] %in% c("Fibroblast") ]
nouse=Reduce(rbind,list( data.frame(den=gis_RGS5,group='Epithelial') ,data.frame(den=gis_RGS51,group='Fibroblast')   ))
pic=list()
pic[[1]]=ggplot(nouse,aes(x=den,fill=group))+geom_density()+mytheme+scale_fill_manual(values = my36colors)+theme(legend.position = 'none')+xlim(-0.4,1.5)+ylim(0,3)

count.one.type.epi$recluster1 <-count.one.type.epi$recluster1 %>%as.factor()
for(i in levels(count.one.type.epi$recluster1)[c(1,2,3)] ){
  gis_RGS53 <- cnv_norm$gis[count.one.type.epi@meta.data[names(cnv_norm$gis),'recluster1']==i]
  nouse=data.frame(den=gis_RGS53,group=i) 
  pic[[i]]=ggplot(nouse,aes(x=den,fill=group))+geom_density()+mytheme+scale_fill_manual(values = my36colors)+xlim(-0.4,1.5)+ylim(0,3)
}
p=CombinePlots(pic,ncol=1)
p

ggsave(paste0(save.pic,'Epi.genomicInsta.pdf'),height = 8,width = 5,p)



save.pic='/data1/huanggy/SNX10_BRCA_project/SC/pic/'

#------DEGS------------------------
Epi$recluster2 <- 'Other'
Epi@meta.data[Epi@meta.data$recluster %in% c('Epi_TREM2'),'recluster2']='Epi_TREM2'
markers=wilcoxauc(Epi,group_by = 'recluster2') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
write.csv(markers,file ='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Epi_TREM2_Meta_T4.csv')



#dotplot--------------------
markers=wilcoxauc(Epi,counts_seu = 'recluster1') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}

DotPlot(counts_seu ,group.by = 'recluster1' ,features =nouse1 ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = 'white', mid = '#fef4d2', high ='purple')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=1, vjust=1,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=60,hjust=1, vjust=1,size = rel(1),color = 'black'))

#Venn------------------------------
Meta_DEG <- read.csv('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Epi_TREM2_Meta_T4.csv',header = T,row.names = 1)
Meta_DEG<-subset(Meta_DEG,subset= group %in%c('Epi_TREM2')) 
Primary_DEG <- read.csv('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Epi_SNX10_Primary.csv',header = T,row.names = 1)
Primary_DEG<-subset(Primary_DEG,subset= group %in%c('Epi_SNX10')) 
CTC_DEG <- read.csv('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Epithelial_Myeloid_CTC.csv',header = T,row.names = 1)
CTC_DEG<-subset(CTC_DEG,subset= group %in%c('Myeloid_CTC')) 

a=intersect(CTC_DEG$feature,Primary_DEG$feature)
intersect(a,Meta_DEG$feature)
Geneset1 <- intersect(a,Meta_DEG$feature)%>%data.frame()
write.csv(Geneset1,file ='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Geneset1.csv')



#-----pathway------
#hallmark
tmp<-Meta_DEG[Meta_DEG$feature%in% Geneset1$. ,]
kegmt<-read.gmt('/data8/huanggy/TCGA_database/data/GSEA_gmt/h.all.v7.0.symbols.gmt') #读gmt文件
kegmt[,1]=capitalize(tolower(str_remove(kegmt[,1],'HALLMARK_')))
hark<-GSEA(tmp$logFC %>% { names(.)=tmp$feature;.} %>% sort(decreasing = T),TERM2GENE = kegmt) #GSEA分析
p=ggplot(hark@result,aes(x=NES,y=reorder(Description,NES)))+geom_segment(aes(yend=Description),xend=0,colour='grey50')+mytheme+geom_point(size=5,aes(colour = -log10(p.adjust)))+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(axis.text.y = element_text(hjust = 1))

#GO
kegmt<-read.gmt('/data8/huanggy/TCGA_database/data/GSEA_gmt/c5.go.v7.4.symbols.gmt') #读gmt文件
kegmt[,1]=kegmt[,1] %>% str_remove('GO_')  %>% tolower()%>% capitalize()
hark<-GSEA(tmp$logFC %>% { names(.)=tmp$feature;.} %>% sort(decreasing = T),TERM2GENE = kegmt) #GSEA分析
p=ggplot(head(hark@result,30),aes(x=NES,y=reorder(Description,NES)))+geom_segment(aes(yend=Description),xend=0,colour='grey50')+mytheme+geom_point(size=5,aes(colour =-log10(p.adjust)))+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(axis.text.y = element_text(hjust = 1))

#KEGG

kegmt<-read.gmt('/data8/huanggy/TCGA_database/data/GSEA_gmt/c2.cp.kegg.v7.4.symbols.gmt') #读gmt文件
kegmt[,1]=kegmt[,1] %>% str_remove('KEGG_')  %>% tolower()%>% capitalize()
hark<-GSEA(tmp$logFC %>% { names(.)=tmp$feature;.} %>% sort(decreasing = T),TERM2GENE = kegmt) #GSEA分析
p=ggplot(head(hark@result,30),aes(x=NES,y=reorder(Description,NES)))+geom_segment(aes(yend=Description),xend=0,colour='grey50')+mytheme+geom_point(size=5,aes(colour = -log10(p.adjust)))+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(axis.text.y = element_text(hjust = 1))

Geneset1 <- intersect(a,Meta_DEG$feature)%>%data.frame()
CELL_LAM <- c('TREM2','LIPA','LPL','CTSB','CTSL','FABP4','FABP5','LGALS1','LGALS3','CD9','CD36')
CR_LAM <- c('TREM2','ACP5','ACP2','LGMN','DAB2','FABP5','LAMP1','LAMP2','LIPA','GPNMB','MSR1','APOE','CD9','FOLR2','CD63','CD68','APOC1')

intersect(CR_LAM,Geneset1$.)
intersect(CELL_LAM,Geneset1$.)
venn_ploy <- venn.diagram(
  x = list(
    Metastasis = Meta_DEG$feature,
    Primary = Primary_DEG$feature,
    CTC = CTC_DEG$feature
  ),
  filename = NULL,
  fill = my36colors[1:3]
)
pdf("/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Geneset1.pdf",width = 5, height = 5) 
grid.draw(venn_ploy)
dev.off()




Epi$recluster2 <- 'Other'
Epi@meta.data[Epi@meta.data$recluster %in% c('Epi_TREM2'),'recluster2']='Epi_TREM2'


markers=wilcoxauc(Epi,group_by = 'recluster2') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}

write.csv(markers,file ='/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Epi_TREM2_Meta_T4.csv')




p1=FeaturePlot(Epi_BrainM,'EPCAM',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p2=FeaturePlot(Epi_BrainM,'KRT19',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p3=FeaturePlot(Epi_BrainM,'CD14',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p4=FeaturePlot(Epi_BrainM,'CD68',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p=plot_grid(p1,p2,p3,p4,ncol = 2, align = "h")
ggsave('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/4_UMAP.pdf',width=7,height=6,p)



p1=FeaturePlot(Epi_BrainM,'EPCAM',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p2=FeaturePlot(Epi_BrainM,'KRT19',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p3=FeaturePlot(Epi_BrainM,'CD14',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p4=FeaturePlot(Epi_BrainM,'CD68',max.cutoff = 1,cols =my36colors)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p=plot_grid(p1,p2,p3,p4,ncol = 2, align = "h")
ggsave('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/4_UMAP.pdf',width=7,height=6,p)






#pie chart--------------

nouse = table(Epi_BrainM$recluster) %>% data.frame()
colnames(nouse) = c('types','num')
nouse$percent = nouse$num/sum(nouse$num)


for (i in seq(nrow(nouse), 1)) {
  if (i == nrow(nouse)) {
    nouse$per.y1[i] = nouse$percent[i]/2 
  }else{
    nouse$per.y1[i] = sum(nouse$percent[(i + 1):nrow(nouse)]) + nouse$percent[i]/2
  }
}


nouse$label = paste(nouse$recluster,'(',round(nouse$percent*100, 2),'%',')', sep = '')


p = ggplot(nouse) +
  geom_bar(aes(x='', percent, fill = label), 
           stat = 'identity', width = .8, color = 'white') +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_manual(values = my36colors) 
ggsave('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/Meta_pie.pdf',width = 7,height =6,p)






#SNX10 expression -------------------------------

p2=FeaturePlot(Epi_BrainM,'SNX10',max.cutoff = 0.5)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p3=FeaturePlot(Epi_BrainM,'CBX3',max.cutoff = 0.5)+mytheme+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red'))
p=plot_grid(p2,p3,ncol = 2, align = "h")
p
ggsave('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/SNX10.pdf',p,height = 4,width=9)

#-----IRGSEA------

CELL_LAM <- c('TREM2','LIPA','LPL','CTSB','CTSL','FABP4','FABP5','LGALS1','LGALS3','CD9','CD36')
CR_LAM <- c('TREM2','ACP5','ACP2','LGMN','DAB2','FABP5','LAMP1','LAMP2','LIPA','GPNMB','MSR1','APOE','CD9','FOLR2','CD63','CD68','APOC1')
kegmt<-read.gmt('/data1/huanggy/GSEA_gmt/c5.go.v7.4.symbols.gmt')
kegmt[,1]=kegmt[,1] %>% str_remove('GO_')  %>% tolower()%>% capitalize()
markers=list()
markers$`Lipid_metabolism`= kegmt[kegmt$term=='Gobp_lipid_metabolic_process',2]#`XXX`用来解决-字符
markers$`CR_LAM`<-CR_LAM
markers$`CELL_LAM`<-CELL_LAM
count.one.type.pathway.Metastasis <- irGSEA.score(object =Epi_BrainM, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                       min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                       species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                       method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                       kcdf = 'Gaussian')
count.one.type.pathway.Primary <- irGSEA.score(object =Onco2021_Epi, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                            min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                            species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                            method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                            kcdf = 'Gaussian')
count.one.type.pathway.CTC <- irGSEA.score(object =counts.seu, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                               min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                               species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                               method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                               kcdf = 'Gaussian')

irGSEA.ridgeplot(object = count.one.type.pathway,method = "AUCell",group.by = 'recluster2',show.geneset = "CR-LAM")+scale_color_manual(values=my36colors)
irGSEA.ridgeplot(object = count.one.type.pathway,method = "AUCell",group.by = 'recluster2',show.geneset = "CELL-LAM")+scale_color_manual(values=my36colors)

irGSEA.ridgeplot(object = count.one.type.pathway,method = "ssgsea",group.by = 'recluster2',show.geneset = "CR-LAM")+scale_color_manual(values=my36colors)
irGSEA.ridgeplot(object = count.one.type.pathway,method = "ssgsea",group.by = 'recluster2',show.geneset = "CELL-LAM")+scale_color_manual(values=my36colors)

p1=irGSEA.halfvlnplot(object = count.one.type.pathway,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "CR-LAM")
p2=irGSEA.halfvlnplot(object = count.one.type.pathway,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "CELL-LAM")
p3=irGSEA.halfvlnplot(object = count.one.type.pathway,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "Lipid-metabolism")
p=plot_grid(p1,p2,p3,ncol =1, align = "h")
p
ggsave('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/LAM_Signature.pdf',p,height = 9,width=5)


counts.seu$recluster1=counts.seu$recluster
counts.seu@meta.data[counts.seu$recluster=='Myeloid' & counts.seu$Tissue=='BloodCTC','recluster1']='Myeloid_CTC'
counts.seu@meta.data[counts.seu$recluster=='Myeloid' & counts.seu$Tissue=='Primary','recluster1']='Myeloid_Primary'
counts.seu@meta.data[counts.seu$recluster%in% c( 'Epithelial_1','Epithelial_2','Epithelial_3','Epithelial_Prolif') & counts.seu$Tissue=='Primary','recluster1']='Epithelial_Primary'
counts.seu@meta.data[counts.seu$recluster%in% c( 'Epithelial_1','Epithelial_2','Epithelial_3','Epithelial_Prolif') & counts.seu$Tissue=='BloodCTC','recluster1']='Epithelial_CTC'


CELL_LAM <- c('TREM2','LIPA','LPL','CTSB','CTSL','FABP4','FABP5','LGALS1','LGALS3','CD9','CD36')
CR_LAM <- c('TREM2','ACP5','ACP2','LGMN','DAB2','FABP5','LAMP1','LAMP2','LIPA','GPNMB','MSR1','APOE','CD9','FOLR2','CD63','CD68','APOC1')
kegmt<-read.gmt('/data4/huanggy/GSEA_gmt/c5.go.v7.4.symbols.gmt')
kegmt[,1]=kegmt[,1] %>% str_remove('GO_')  %>% tolower()%>% capitalize()
markers=list()
markers$Lipid= kegmt[kegmt$term=='Gobp_lipid_metabolic_process',2]
markers$CR_LAM<-CR_LAM
markers$CELL_LAM<-CELL_LAM
count.one.type.pathway.Metastasis <- irGSEA.score(object =Epi_BrainM, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                                  min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                                  species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                                  method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                                  kcdf = 'Gaussian')
count.one.type.pathway.Primary <- irGSEA.score(object =Onco2021_Epi, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                               min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                               species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                               method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                               kcdf = 'Gaussian')
count.one.type.pathway.CTC <- irGSEA.score(object =counts.seu, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                           min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                           species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                           method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                           kcdf = 'Gaussian')
p=DimPlot(counts.seu,group.by = 'recluster1',cols = my36colors)+mytheme
ggsave('/data4/huanggy/CTC_recluster1.pdf',p,height = 5,width=5)

p1=irGSEA.halfvlnplot(object = count.one.type.pathway.CTC,
                      method = "AUCell",group.by = 'recluster1',
                      show.geneset = "CR-LAM")+mytheme
p2=irGSEA.halfvlnplot(object = count.one.type.pathway.Primary,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "CR-LAM")+mytheme
p3=irGSEA.halfvlnplot(object = count.one.type.pathway.Metastasis,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "CR-LAM")+mytheme

p4=irGSEA.halfvlnplot(object = count.one.type.pathway.CTC,
                      method = "AUCell",group.by = 'recluster1',
                      show.geneset = "CELL-LAM")+mytheme
p5=irGSEA.halfvlnplot(object = count.one.type.pathway.Primary,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "CELL-LAM")+mytheme
p6=irGSEA.halfvlnplot(object = count.one.type.pathway.Metastasis,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "CELL-LAM")+mytheme




p7=irGSEA.halfvlnplot(object = count.one.type.pathway.CTC,
                      method = "AUCell",group.by = 'recluster1',
                      show.geneset = "Lipid")+mytheme
p8=irGSEA.halfvlnplot(object = count.one.type.pathway.Primary,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "Lipid")+mytheme
p9=irGSEA.halfvlnplot(object = count.one.type.pathway.Metastasis,
                      method = "AUCell",group.by = 'recluster',
                      show.geneset = "Lipid")+mytheme

p=plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol =3, align = "h")
p

ggsave('/data4/huanggy/SNX10_BRCA_project/SNX10_T4_SC/pic/LAM_Signature.pdf',p,height = 9,width=5)


p=plot_grid(p1,p4,p7,ncol =1, align = "h")
p

ggsave('/data4/huanggy/CTC.pdf',p,height = 9,width=7)
p=plot_grid(p2,p5,p8,ncol =1, align = "h")
p

ggsave('/data4/huanggy/Primary.pdf',p,height = 9,width=11)


p=plot_grid(p3,p6,p9,ncol =1, align = "h")
p

ggsave('/data4/huanggy/BrainM.pdf',p,height = 9,width=4)


#TCGA_COR-------------------

nouse <- data.frame(cancertype=c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SRAC','SKCM','STAD','TCGT','THCA','THYM','UCEC','UCS','UVM')
                    ,cor=c(0.75,0.34,0.54,0.36,0.64,0.16,0.21,0.44,0.23,0.42,0.13,0.73,0.41,0.19,0.13,0.59,0.23,0.3,0.094,0.39,0.28,0.4,0.58,0.2,0.28,0.66,0.62,0.62,0.26,0.43,0.37,0.35,0.78)
                    ,p=c(0.001,0.001,0.001,0.001,0.001,0.008,0.16,0.001,0.004,0.001,0.29,0.001,0.001,0.011,0.003,0.001,0.001,0.001,0.39,0.001,0.001,0.001,0.001,0.061,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.008,0.001))
p=ggplot(nouse,aes(x=reorder(cancertype,cor),y=cor,color=-log10(p)))+geom_point(size=7)+scale_color_gradientn(colors = c('#cacaca','#ffbf8750','red')) +mytheme
ggsave('/data1/huanggy/SNX10_BRCA_project/33cancertype_cor_CBX3_SNX10.pdf',p,height = 4,width=12)

#Geneset1 enrichGo--------------
gene.name<-bitr(Geneset1$., fromType = "SYMBOL", 
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db,drop =  FALSE )
gene.name<-gene.name[!is.na(gene.name$ENSEMBL),]
gene.name<-merge(x=tmp,y=gene.name,by.x='name',by.y='SYMBOL')
geneList<-gene.name$logFC
names(geneList) = gene.name$ENTREZID
geneList= sort(geneList,decreasing = T)
eG <- enrichGO(gene = gene.name$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.01, 
               qvalueCutoff = 0.01, 
               ont="all", 
               readable =T)
write.csv(eG,'/data5/huanggy/Geneset1_enrichgo.csv')
eG <- read.csv('/data1/huanggy/SNX10_BRCA_project/Geneset1_enrichgo.csv',row.names = 1)
select <-c('GO:0002757',
           'GO:0006909',
           'GO:0070661',
           'GO:0097529',
           'GO:0030041',
           'GO:0022409',
           'GO:0019216',
           'GO:0006801',
           'GO:0045055',
           'GO:0043112',
           'GO:0007015',
           'GO:0042116',
           'GO:0043410',
           'GO:0002253',
           'GO:0060326'
)

eG <- separate(data=eG, col=GeneRatio,into = c("GR1", "GR2"), sep = "/")
eG <- separate(data=eG, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
eG <- mutate(eG, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
p <- ggplot(eG[select,],aes(enrichment_factor,reorder(Description,enrichment_factor))) + 
  geom_point(aes(size=Count,color=-log10(pvalue))) +
  scale_color_gradientn(colors = c('royalblue1','red')) + 
  labs(color="-log10(pvalue)",size="Count",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + mytheme

nouse <- eG[select,]
p <- ggplot(eG[select,],aes(-log10(p.adjust),reorder(Description,-log10(p.adjust)),fill=-log10(p.adjust))) + geom_bar(stat='identity')+scale_fill_gradientn(colors = c('royalblue1','red'))  + mytheme

ggsave('/data1/huanggy/SNX10_BRCA_project/Geneset_enrichgo.pdf',height = 6,width = 8,p)


#CTC_SC--------------------
markers=wilcoxauc(counts.seu,group_by = 'recluster1') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}

p=DotPlot(counts.seu ,group.by = 'recluster1' ,features =c('CD3D','CD3E','CD14','LYZ','LIPA','TYROBP','COL1A1','COL6A1','EPCAM','KRT8','KRT19','KET8') ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = 'white', mid = '#fef4d2', high ='red')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=1, vjust=1,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=60,hjust=1, vjust=1,size = rel(1),color = 'black'))

ggsave('/data4/huanggy/Fusion/CTC_SC_Dot.pdf',width=6,height=4,p)

#find SNX10------
Geneset1 <- read.csv('/data1/huanggy/SNX10_BRCA_project/SNX10_T4_SC/data/Geneset1.csv')
primary <- read.csv('/data1/huanggy/SNX10_BRCA_project/Primary_Epithelial_Myeloid_CTC.csv',row.names = 1)
primary <- subset(primary,subset=group=='Myeloid_CTC')

CTC <- read.csv('/data1/huanggy/SNX10_BRCA_project/Epithelial_Myeloid_CTC.csv',row.names = 1)
CTC <- subset(CTC,subset=group=='Myeloid_CTC')

resRNA_1 <- primary
resRNA_2 <- CTC
resRNA_1_change <- resRNA_1[resRNA_1$feature %in% intersect(resRNA_1$feature,resRNA_2$feature),]
resRNA_2_change <- resRNA_2[resRNA_2$feature%in% intersect(resRNA_1$feature,resRNA_2$feature),]
resRNA_2_change<- resRNA_2_change[rownames(resRNA_1_change),]

tmp.pic=data.frame(primary=resRNA_1_change$logFC,CTC=resRNA_2_change$logFC,feature=resRNA_1_change$feature)
tmp.pic[,'group']='No.sig'
tmp.pic[tmp.pic$primary> 0.58 & tmp.pic$CTC> 0.58,'group']='Myeloid_CTC.Both.up'
tmp.pic[tmp.pic$primary< -0.58 & tmp.pic$CTC< -0.58,'group']='Epithelial.Both.up'
tmp.pic[tmp.pic$primary> 0.58 & tmp.pic$CTC< -0.58,'group']='primary.Up'
tmp.pic[tmp.pic$primary< -0.58 & tmp.pic$CTC> 0.58,'group']='CTC.Up'

tophit= tmp.pic[tmp.pic$group!='No.sig',]
CELL_LAM <- c('TREM2','LIPA','LPL','CTSB','CTSL','FABP4','FABP5','LGALS1','LGALS3','CD9','CD36')

CR_LAM <- c('TREM2','ACP5','ACP2','LGMN','DAB2','FABP5','LAMP1','LAMP2','LIPA','GPNMB','MSR1','APOE','CD9','FOLR2','CD63','CD68','APOC1')

tophit1 <- tmp.pic[tmp.pic$feature %in%intersect(tmp.pic$feature,c('SNX10',CR_LAM,CELL_LAM)),]
p=ggplot(tmp.pic,aes(x=primary,y=CTC,color=group))+geom_point(size=5)+mytheme+scale_color_manual(values = c(my36colors[1],my36colors[2],my36colors[3],my36colors[4],'#66666640'))+
  geom_vline(aes(xintercept=-0.58),colour="darkgrey", linetype="dashed")+
  geom_vline(aes(xintercept=0.58),colour="darkgrey", linetype="dashed")+
  geom_hline(aes(yintercept=0.58),colour="darkgrey", linetype="dashed")+
  geom_hline(aes(yintercept=-0.58),colour="darkgrey", linetype="dashed")+
  geom_text_repel(data = tophit1,aes(label = feature),size = 3,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"),segment.color='black',colour = "black")
ggsave('/data1/huanggy/SNX10_BRCA_project/Myeloid_CTC_Epithelial_tmppic.pdf',width =8,height =6.5,p)
write.csv(tmp.pic,'/data1/huanggy/SNX10_BRCA_project/Myeloid_CTC_Epithelial_tmppic.csv')



#Dimplot------------------
counts$recluster2 <- counts$recluster
counts@meta.data[counts@meta.data$recluster1 %in% c('Epi_SNX10'),'recluster2'] <- 'Fusioncell'
p=DimPlot(counts,group.by = 'recluster2',cols  = c('#ffbf8750','#ffbf8750','red','#ffbf8750','yellow','#ffbf8750','royalblue','#ffbf8750','#ffbf8750'))+mytheme
ggsave('/data4/huanggy/Fusion/Primary_recluster2.pdf',height = 5,width=6,p)
p=DimPlot(counts,group.by = 'recluster2',cols  = my36colors)+mytheme

ggsave('/data4/huanggy/Fusion/Primary_recluster.pdf',height = 5,width=6,p)

Epithelial_primary
p=DimPlot(Epithelial_primary,group.by = 'recluster',cols  = my36colors)+mytheme
ggsave('/data4/huanggy/Fusion/Primary_recluster1.pdf',height = 5,width=6,p)

counts.seu$group <- 'NonMutation'
counts.seu@meta.data[c('B2_BloodCTC_Bre_X_3','B2_BloodCTC_Bre_X_6','B2_BloodCTC_Bre_X_7','B2_BloodCTC_Bre_X_8','B2_BloodCTC_Bre_X_9','B2_BloodCTC_Bre_X_11','B2_BloodCTC_Bre_X_13','B2_BloodCTC_Bre_X_14'),'group'] <- 'Mutation'


p=DimPlot(counts.seu,group.by = 'recluster1',cols  = my36colors)+mytheme

ggsave('/data4/huanggy/Fusion/CTC_recluster1.pdf',height = 5,width=6,p)


p=DimPlot(counts.seu,group.by = 'group',cols  = c('red3','#cacaca'))+mytheme

ggsave('/data4/huanggy/Fusion/CTC_Mye_Mut.pdf',height = 5,width=6,p)

p=DimPlot(counts.Meta,group.by = 'recluster',cols  = my36colors)+mytheme
counts.Meta$recluster1 <- counts.Meta$recluster
counts.Meta@meta.data[rownames(Epi_BrainM@meta.data[Epi_BrainM$recluster %in% c('Epi_TREM2'),]),'recluster1'] <-'Fusioncell'
counts.Meta@meta.data[rownames(Epi_Meta@meta.data[Epi_Meta$recluster %in% c('Epi_TREM2'),]),'recluster1'] <-'Fusioncell'

p=DimPlot(counts.Meta,group.by = 'recluster1',cols  = my36colors)+mytheme
ggsave('/data4/huanggy/Fusion/Meta_recluster.pdf',height = 5,width=6,p)
p=DimPlot(counts.Meta,group.by = 'recluster',cols  = my36colors)+mytheme
ggsave('/data4/huanggy/Fusion/Meta_recluster1.pdf',height = 5,width=6,p)

p=DimPlot(counts.Meta,group.by = 'recluster1',cols  = c('#ffbf8750','red','#ffbf8750','yellow','royalblue','#ffbf8750','#ffbf8750'))+mytheme
ggsave('/data4/huanggy/Fusion/Meta_recluster2.pdf',height = 5,width=6,p)

gene.choose1<-c('KRT8','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','CD3D','CD3E')

DotPlot(counts.seu,group.by = 'recluster1' ,features = gene.choose1  ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = 'white', mid = '#fef4d2', high ='red')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=1, vjust=1,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=1, vjust=1,size = rel(1),color = 'black'))+theme(legend.position = 'none')



#density---------------------
p1=ggplot(counts.seu@meta.data, aes(x = nCount_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p2=ggplot(counts.seu@meta.data, aes(x = nFeature_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p=plot_grid(p1,p2,ncol=1,align='h')
ggsave('/data4/huanggy/Fusion/quality_CTC.pdf',height = 5,width=5,p)



p1=ggplot(counts@meta.data, aes(x = nCount_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p2=ggplot(counts@meta.data, aes(x = nFeature_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p=plot_grid(p1,p2,ncol=1,align='h')
ggsave('/data4/huanggy/Fusion/quality_primary.pdf',height = 5,width=5,p)


p1=ggplot(Epi_BrainM@meta.data, aes(x = nCount_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p2=ggplot(Epi_BrainM@meta.data, aes(x = nFeature_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p=plot_grid(p1,p2,ncol=1,align='h')
ggsave('/data4/huanggy/Fusion/quality_CTC.pdf',height = 5,width=3,p)



p1=ggplot(counts.Meta@meta.data, aes(x = nCount_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p2=ggplot(counts.Meta@meta.data, aes(x = nFeature_RNA,fill=my36colors[1])) +geom_density(alpha=.3) +scale_fill_manual(values = my36colors)+mytheme
p=plot_grid(p1,p2,ncol=1,align='h')
ggsave('/data4/huanggy/Fusion/quality_Meta.pdf',height = 5,width=5,p)

#Violin Diagram--------------

Gene.choose<-c('CD79B','MS4A1','PECAM1','VWF','KRT8','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','CPA3','MS4A2','TNFRSF17','MZB1','CD3D','CD3E')


p=StackedVlnPlot(counts, Gene.choose, 'recluster2',color=my36colors,pt.size=0)
ggsave('/data4/huanggy/Fusion/Primary_violin.pdf',height = 8,width=6,p)


Epithelial_primary
Gene.choose<-c('CD79B','MS4A1','PECAM1','VWF','KRT8','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','CPA3','MS4A2','TNFRSF17','MZB1','CD3D','CD3E')


p=StackedVlnPlot(Epithelial_primary, Gene.choose, 'recluster',color=my36colors,pt.size=0)
ggsave('/data4/huanggy/Fusion/Primary_violin1.pdf',height = 8,width=6,p)


Gene.choose<-c('KRT8','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','CD3D','CD3E')


p=StackedVlnPlot(counts.seu, Gene.choose, 'recluster1',color=my36colors,pt.size=0)
ggsave('/data4/huanggy/Fusion/Primary_violin.pdf',height = 8,width=6,p)



Gene.choose<-c('PECAM1','VWF','KRT19','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','OLIG1','MBP','CD3D','CD3E')


p=StackedVlnPlot(counts.Meta, Gene.choose, 'recluster1',color=my36colors,pt.size=0)
ggsave('/data4/huanggy/Fusion/Meta_violin.pdf',height = 8,width=6,p)



Gene.choose<-c('PECAM1','VWF','KRT19','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','OLIG1','MBP','CD3D','CD3E')


p=StackedVlnPlot(counts.Meta, Gene.choose, 'recluster',color=my36colors,pt.size=0)
ggsave('/data4/huanggy/Fusion/Meta_violin1.pdf',height = 8,width=6,p)


