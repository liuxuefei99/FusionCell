.libPaths('/data1/liuxf/R_package4')
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
library(SingleCellExperiment)
library(smoother)
library(plotrix)
library(ggpubr)
library(RColorBrewer)
library(Seurat)

library(magrittr)
library(stringr)
library(ggalluvial)
library(ggforce)
library(org.Mm.eg.db)
library(AUCell)
library(ape)
library(Hmisc)
library(survival)
library(survminer)
library(survivalROC)
library(ggimage)
library(gridGraphics)
library(cowplot)


my36colors <-c('#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1",'red','blue','green','red','red')
my36colors1<-c('#594C57',"#8F6D5F", '#EAE5E3','#DB9C5E',"#420B2F","#C25160","#CC5D20",'#2E59A7','#80A492','#32788A','#E1D279','#A64036')

my36colors2 <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"))  


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

main.path='/data2/chengyx/cellfusion'
save.data='/data2/chengyx/cellfusion/data/'
save.pic='/data2/chengyx/cellfusion/pic/'

#===+=====载入函数=====+=====


#-----loading -------
set.seed(123)

counts.breast=readRDS(paste0(save.data,'BRCA.count.rds'))
counts.breast=counts.breast[,colnames(counts.breast) %>% str_count('Primary') %>% as.logical()]
counts.CTC=readRDS(paste0(save.data,'BrianMeta.count.rds'))
counts.CTC=counts.CTC[,colnames(counts.CTC) %>% str_count('Bre') %>% as.logical()]

counts.breast=cbind(counts.breast,counts.CTC)

tmp.file=read.delim(paste0(save.data,'GSE111065.csv'),sep=',',header=T,row.names = 1)
counts.CTC1=tmp.file[,colnames(tmp.file) %>% str_count('_BloodCTC') %>% as.logical()]
colnames(counts.CTC1) = colnames(counts.CTC1)  %>% str_remove('BRCA_')

counts.breast=cbind(counts.breast,counts.CTC1)

#tmp.file=read.delim(paste0(save.data,'GSE109761.csv'),sep=',',header=T,row.names = 1)

saveRDS(counts.breast,paste0(save.data,'BrianBreast_Breast.count.rds'))

counts.seu <- CreateSeuratObject(counts = counts.breast)
counts.seu[["percent.mt"]] = PercentageFeatureSet(counts.seu, pattern = "^[Mm]T-")
counts.seu[["RP.mt"]] = PercentageFeatureSet(counts.seu, pattern = "^RP[LS]\\w")
counts.seu <- subset(counts.seu, subset = (nFeature_RNA > 200  & percent.mt < 40 ))


counts.seu <- NormalizeData(object =counts.seu, normalization.method = "LogNormalize", scale.factor = 10000)
counts.seu <- FindVariableFeatures(object = counts.seu, selection.method = "vst", nfeatures = 2000)
counts.seu <- ScaleData(counts.seu, vars.to.regress = c("nCount_RNA"))
counts.seu <- RunPCA(object =counts.seu, features = VariableFeatures(object = counts.seu))

message=rownames(counts.seu@meta.data) %>% str_split_fixed('_',n=7) %>% data.frame()

counts.seu@meta.data[,'Patient']=message[,1]
counts.seu@meta.data[,'Tissue']=message[,2]
counts.seu@meta.data[,'Treat']=message[,4]

counts.seu <- RunHarmony(counts.seu, group.by.vars = c("Patient"), dims.use = 1:33, verbose = T)
counts.seu <- FindNeighbors(object = counts.seu, dims = 1:35, reduction = "harmony")
counts.seu <- FindClusters(object = counts.seu, resolution = 0.5)
counts.seu <- RunUMAP(object = counts.seu, dims = 1:30, reduction = "harmony")

name=c('0'='Epithelial_1','1'='T cell','2'='Myeloid','3'='Epithelial_1','4'='Epithelial_1','5'='Epithelial_2','6'='Epithelial_Prolif','8'='Fibroblast','9'='Epithelial_3')
counts.seu$recluster=name[counts.seu$seurat_clusters %>% as.character()]

saveRDS(counts.seu,paste0(save.data,'counts.breast2Brian.rds'))

p=DimPlot(counts.seu,group.by = 'recluster',pt.size = 2)+mytheme+scale_color_manual(values = my36colors2)
ggsave(paste0(save.pic,'BRCA.CTC.cluster.pdf'),width = 6,height = 5 )

Gene.choose=c('KRT8','EPCAM','COL1A1','COL6A1','CD14','FCGR3A','CD3D','CD3E')
p=StackedVlnPlot(counts.seu, Gene.choose, 'recluster',color=my36colors2,pt.size=0)
ggsave(paste0(save.pic,'BRCA.CTC.markers.pdf'),width = 4,height = 5 )



#------tissue-------
p=DimPlot(counts.seu,group.by = 'Tissue',pt.size = 2)+mytheme+scale_color_manual(values = c('#11111102','red',rep('#11111102',9)))
p=DimPlot(counts.seu,group.by = 'Tissue',pt.size = 2)+mytheme+scale_color_manual(values = c('red',rep('#11111110',9)))
ggsave(paste0(save.pic,'BRCA.CTC.UMAP.pdf'),width = 6,height = 5 )

#--------markers--------
counts.seu$recluster1=counts.seu$recluster
counts.seu@meta.data[counts.seu$recluster=='Myeloid' & counts.seu$Tissue=='BloodCTC','recluster1']='Myeloid_CTC'
counts.seu@meta.data[counts.seu$recluster=='Myeloid' & counts.seu$Tissue=='Primary','recluster1']='Myeloid_Primary'
counts.seu@meta.data[counts.seu$recluster%in% c( 'Epithelial_1','Epithelial_2','Epithelial_3','Epithelial_Prolif') & counts.seu$Tissue=='Primary','recluster1']='Epithelial_Primary'
counts.seu@meta.data[counts.seu$recluster%in% c( 'Epithelial_1','Epithelial_2','Epithelial_3','Epithelial_Prolif') & counts.seu$Tissue=='BloodCTC','recluster1']='Epithelial_CTC'


p=DimPlot(counts.seu,group.by = 'recluster1',pt.size = 2)+mytheme+scale_color_manual(values = my36colors2)
ggsave(paste0(save.pic,'BRCA.CTC.cluster1.pdf'),width = 6,height = 5 )


#-----compare epiCTC and myeloidCTC------

C1.barcode=rownames(counts.seu@meta.data[counts.seu$recluster1=='Myeloid_CTC',])
C2.barcode=rownames(counts.seu@meta.data[counts.seu$recluster1=='Epithelial_CTC',])
now_meta <- data.frame(row.names=c(C1.barcode,C2.barcode), condition=counts.seu@meta.data[c(C1.barcode,C2.barcode),'recluster1'] )
y <- DGEList(counts=as.matrix(counts.seu@assays$RNA@counts[,rownames(now_meta)]),group=now_meta$condition )
y <- y[rowSums(edgeR::cpm(y)>1) >= 2, , keep.lib.sizes=FALSE] %>% calcNormFactors() %>% estimateDisp( model.matrix(~now_meta$condition))
et <- exactTest(y)
res <- et$table %>% dplyr::mutate( name=rownames(.) ) %>% dplyr::mutate(logpvalue= -log10(.$PValue) ) %>%  { .[! .$name %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]} %>% { .[is.infinite(.$logpvalue),'logpvalue'] = max(.$logpvalue[.$logpvalue!=max(.$logpvalue)]) ;. } %>%  {.[.$logpvalue>500,'logpvalue']=500 ;.}
res = res %>% {.[.$logFC>1.5&.$logpvalue>2,'color']='Myeloid_CTC';.}%>% {.[.$logFC< -1.5&.$logpvalue>2,'color']='Epithelial_CTC';. } %>% {.[is.na(.$color) ,'color']= 'No Sig';.}
tophit =res[res$name %in% c('TREM2','CCR5','CD163','FABP4','FCGR1A','MMP9','MARCO','S100A8','CD86','FCN1','LIPA','LPL','CTSB','FABP5','CD36','LGALS3','LYZ','SNX10'),]


xmax=max(abs(min(res$logFC)),max(res$logFC))

p <- ggplot(res, aes(x = logFC, y = logpvalue)) +xlim(-xmax,xmax)+geom_point(aes(color = color),size=2)+
  scale_color_manual(values = c("#E41A1C60","#377EB860","grey")) +mytheme+
  geom_text_repel(data = tophit,aes(label = name),size = 3,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"),segment.color='black',colour = "black")+ 
  geom_vline(aes(xintercept=-1.5),colour="darkgrey", linetype="dashed")+geom_vline(aes(xintercept=1.5),colour="darkgrey", linetype="dashed") +geom_hline(aes(yintercept=2),colour="darkgrey", linetype="dashed")

ggsave(paste0(save.pic,'BRCA.CTC-2DEGS.pdf'),width =8,height =8,p)



#----Infercnv-----

#-------infercnv：Reference：Normal Epi--------

Epi_count.cnv=counts.seu
cell.keep=c(rownames(counts.seu@meta.data[counts.seu$Tissue=='BloodCTC',]) , rownames(counts.seu@meta.data[counts.seu$recluster %in% c('Myeloid','Tcell') & counts.seu$Tissue!='BloodCTC',]))
Epi_count.cnv=Epi_count.cnv %>% subset(cells=cell.keep)


now.meta=Epi_count.cnv@meta.data
system(paste0('mkdir -p ',save.data,'/data/infercnv/BloodCTC_Ref_NormalMyeloid_BRCA'))
Epi_count.cnv$recluster1=Epi_count.cnv$recluster1 %>% as.character()
Epi_count.cnv@meta.data[Epi_count.cnv$recluster%in% c('Myeloid','Tcell') & Epi_count.cnv$Tissue!='BloodCTC','recluster1']='Normal'



ann=read.csv(paste0(save.data,'/now_gene_pos.txt'),header=F,sep='\t')
count=Epi_count.cnv@assays$RNA@counts
intersectgene=intersect(rownames(count),ann$V1)
ann_new=ann[ann$V1 %in% intersectgene,]
count_new=count[ann_new$V1,] %>% data.frame()

meta=data.frame(row.names = rownames(Epi_count.cnv@meta.data),recluster=Epi_count.cnv$recluster1)

ann_new=data.frame(row.names = ann_new$V1,ann_new[2:4])

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_new,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=ann_new,
                                    ref_group_names=c("Normal"),
                                    chr_exclude=c('chrM')) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0(save.data,'/data/infercnv/BloodCTC_Ref_NormalMyeloid_BRCA'),  # 输出文件夹
                             cluster_by_groups=T,# 聚类
                             num_threads=10,
                             noise_filter=0.2,
                             denoise=T, #去噪
                             HMM=T,no_plot = F) # 是否基于HMM预测CNV

add_to_seurat(seurat_obj = Epi_count.cnv, paste0(save.data,'/data/infercnv/BloodCTC_Ref_NormalMyeloid_BRCA/'), top_n = 10,bp_tolerance = 2e+06)


