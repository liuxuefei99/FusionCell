


Fusion_mut <-read.table('/data4/huanggy/WES_1/3.mutation/Fusion/Fusion.hg38_multianno.txt', sep = "\t", header = T,encoding = "utf-8")

LAM_mut <-read.table('/data4/huanggy/WES_1/3.mutation/LAM/LAM.hg38_multianno.txt', sep = "\t", header = T)
WT_mut <-read.table('/data4/huanggy/WES_1/3.mutation/NLS/NLS.hg38_multianno.txt', sep = "\t", header = T)



Fusion_mut = annovarToMaf(annovar = "/data4/huanggy/WES_1/3.mutation/Fusion/Fusion.hg38_multianno.txt", 
                               Center = 'NA', 
                               refBuild = 'hg38', 
                               tsbCol = 'Tumor_Sample_Barcode', 
                               table = 'refGene',
                               sep = "\t")
Fusion_mut <- Fusion_mut[,c(1,2,3,12,15,127)]
colnames(Fusion_mut)[6]<- 'DP'
Fusion_mut$snp <- paste0(Fusion_mut$Start_Position,'_',Fusion_mut$Gene.refGene)
Fusion_mut$DP <- Fusion_mut$DP %>%as.numeric()


LAM_mut = annovarToMaf(annovar = "/data4/huanggy/WES_1/3.mutation/LAM/LAM.hg38_multianno.txt", 
                          Center = 'NA', 
                          refBuild = 'hg38', 
                          tsbCol = 'Tumor_Sample_Barcode', 
                          table = 'refGene',
                          sep = "\t")
LAM_mut <- LAM_mut[,c(1,2,3,12,15,127)]
colnames(LAM_mut)[6]<- 'DP'
LAM_mut$snp <- paste0(LAM_mut$Start_Position,'_',LAM_mut$Gene.refGene)
LAM_mut$DP <- LAM_mut$DP %>%as.numeric()



WT_mut = annovarToMaf(annovar = "/data4/huanggy/WES_1/3.mutation/NLS/NLS.hg38_multianno.txt", 
                          Center = 'NA', 
                          refBuild = 'hg38', 
                          tsbCol = 'Tumor_Sample_Barcode', 
                          table = 'refGene',
                          sep = "\t")
WT_mut <- WT_mut[,c(1,2,3,12,15,127)]
colnames(WT_mut)[6]<- 'DP'
WT_mut$snp <- paste0(WT_mut$Start_Position,'_',WT_mut$Gene.refGene)
WT_mut$DP <- WT_mut$DP %>%as.numeric()









WT_mut1 <- WT_mut
LAM_mut1 <- LAM_mut
Fusion_mut1 <- Fusion_mut

setdiff(intersect(Fusion_mut1$snp,LAM_mut1$snp),WT_mut1$snp)

all3 <- intersect(intersect(Fusion_mut$snp,LAM_mut$snp),WT_mut$snp)
Fusion_LAM <- intersect(Fusion_mut$snp,LAM_mut$snp)
Fusion_WT <- intersect(Fusion_mut$snp,WT_mut$snp)

unique_Fusion_LAM <- setdiff(Fusion_LAM,all3)
unique_Fusion_WT <- setdiff(Fusion_WT,all3)


cnv_LAM <- Fusion_mut[Fusion_mut$snp %in%unique_Fusion_LAM,]
write.csv(cnv_LAM,'/data4/huanggy/WES_1/cnv_LAM.csv')
cnv_WT <- Fusion_mut[Fusion_mut$snp %in%unique_Fusion_WT,]
write.csv(cnv_WT,'/data4/huanggy/WES_1/cnv_WT.csv')


setdiff(cnv_LAM$snp,cnv_WT$snp)


venn_ploy <- venn.diagram(
  x = list(
    Fusion = Fusion_mut$snp,
    LAM = LAM_mut$snp,
    WT = WT_mut$snp
  ),
  filename = NULL,
  fill = my36colors[1:3]
)


pdf("/data4/huanggy/Fusion/cnv_event.pdf",width = 5, height = 5) 
grid.draw(venn_ploy)
dev.off()



DriverGenes <- read_tsv('/data4/huanggy/WES_1/IntOGen-DriverGenes.tsv')

intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol) 
intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) 


WT_unique <- intersect(str_split_fixed(setdiff(Fusion_WT,Fusion_LAM),'_',2)[,2],DriverGenes$Symbol) 
WT_unique_site <- Fusion_mut[Fusion_mut$snp %in% setdiff(Fusion_WT,Fusion_LAM),]
LAM_unique <- intersect(str_split_fixed(setdiff(Fusion_LAM,Fusion_WT),'_',2)[,2],DriverGenes$Symbol) 
LAM_unique_site <- Fusion_mut[Fusion_mut$snp %in% setdiff(Fusion_LAM,Fusion_WT),]

LAM_unique_site1 <- LAM_unique_site[LAM_unique_site$DP > 100,]

intersect(LAM_unique_site1$Gene.refGene,DriverGenes$Symbol) 

intersect(str_split_fixed(Fusion_WT,'_',2)[,2],DriverGenes$Symbol) 

Fusion_LAM <- intersect(Fusion_mut$snp,LAM_mut$snp)
Fusion_WT <- intersect(Fusion_mut$snp,WT_mut$snp)
setdiff(Fusion_LAM,Fusion_WT)



setdiff(intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol),intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) )


setdiff(intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) ,intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol))

str_split_fixed(setdiff(cnv_LAM$snp,cnv_WT$snp),'_',2)[,2]  

intersect(str_split_fixed(setdiff(cnv_LAM$snp,cnv_WT$snp),'_',2)[,2],DriverGenes$Symbol)
setdiff(intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) ,intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol))





#部分基因
WT_mut1 <- WT_mut[WT_mut$DP > 400,]
LAM_mut1 <- LAM_mut[LAM_mut$DP > 100,]
Fusion_mut1 <- Fusion_mut[Fusion_mut$DP >400,]

setdiff(intersect(Fusion_mut1$snp,LAM_mut1$snp),WT_mut1$snp)

all3 <- intersect(intersect(Fusion_mut1$snp,LAM_mut1$snp),WT_mut1$snp)
Fusion_LAM <- intersect(Fusion_mut1$snp,LAM_mut1$snp)
Fusion_WT <- intersect(Fusion_mut1$snp,WT_mut1$snp)

unique_Fusion_LAM <- setdiff(Fusion_LAM,all3)
unique_Fusion_WT <- setdiff(Fusion_WT,all3)


cnv_LAM <- Fusion_mut[Fusion_mut1$snp %in%unique_Fusion_LAM,]
write.csv(cnv_LAM,'/data4/huanggy/WES_1/cnv_LAM.csv')
cnv_WT <- Fusion_mut[Fusion_mut1$snp %in%unique_Fusion_WT,]
write.csv(cnv_WT,'/data4/huanggy/WES_1/cnv_WT.csv')


setdiff(cnv_LAM$snp,cnv_WT$snp)


venn_ploy <- venn.diagram(
  x = list(
    Fusion = Fusion_mut$snp,
    LAM = LAM_mut$snp,
    WT = WT_mut$snp
  ),
  filename = NULL,
  fill = my36colors[1:3]
)


pdf("/data4/huanggy/Fusion/cnv_event.pdf",width = 5, height = 5) 
grid.draw(venn_ploy)
dev.off()


library(tidyverse)
DriverGenes <- read_tsv('/data4/huanggy/WES_1/IntOGen-DriverGenes.tsv')

intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol) 
intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) 


intersect(str_split_fixed(Fusion_LAM,'_',2)[,2],DriverGenes$Symbol) 
intersect(str_split_fixed(Fusion_WT,'_',2)[,2],DriverGenes$Symbol) 

intersect(str_split_fixed(unique_Fusion_LAM,'_',2)[,2],DriverGenes$Symbol) 
intersect(str_split_fixed(unique_Fusion_WT,'_',2)[,2],DriverGenes$Symbol) 



setdiff(intersect(Fusion_mut1$snp,LAM_mut1$snp),WT_mut1$snp)





setdiff(intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol),intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) )


setdiff(intersect(cnv_LAM$Gene.refGene,DriverGenes$Symbol) ,intersect(cnv_WT$Gene.refGene,DriverGenes$Symbol))

str_split_fixed(setdiff(cnv_LAM$snp,cnv_WT$snp),'_',2)[,2]  

intersect(str_split_fixed(setdiff(cnv_LAM$snp,cnv_WT$snp),'_',2)[,2],DriverGenes$Symbol)
