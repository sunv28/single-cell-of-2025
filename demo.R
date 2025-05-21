#### loading package ####
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(reshape2)
library(patchwork)

####loading DATA#####
scRNA.data <- Read10X(data.dir = "./demo_data_1/filtered_gene_bc_matrices/hg19/")

scRNA <- CreateSeuratObject(counts = scRNA.data, 
                            project = "pbmc3k",  
                            min.cells = 3,               # Include features where at least this many cells are detected 
                            min.features = 200)          # Include cells where at least this many features are detected
scRNA

####QC metrics and filter####

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

# calculate erythrocyte QC metrics
## HB.genes(including human,mouse,Rat)
HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ","Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-y", "Hbb-bh1", "Hbe1", "Hbg1", "Hbg2", "Hbz","Hba1", "Hba2", "Hbb", "Hbd", "Hbe1", "Hbg1", "Hbg2", "Hbq1", "Hbz")


## Match the expression matrix and select the HB.genes
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 


## calculate erythrocyte QC metrics
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 


scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




#NormalizeData######

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)



####VariableFeatures#####
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)



top10 <- head(VariableFeatures(scRNA), 10)

##scale data####

all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)



###pca#####
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
print(scRNA[["pca"]], dims = 1:5, nfeatures = 5)
##cell cycle#####
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes



RidgePlot(scRNA, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

scRNA <- RunPCA(scRNA, features = c(s.genes, g2m.genes))
DimPlot(scRNA)

?AddModuleScore()

scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))


scRNA <- RunPCA(scRNA, features = c(s.genes, g2m.genes))
DimPlot(scRNA)

####SCTransform#####
scRNA <- SCTransform(scRNA, vars.to.regress = "percent.mt", verbose = FALSE)

a <- scRNA@assays[["RNA"]]@layers[["scale.data"]]

a <- scRNA[["RNA"]]$scale.data
colnames(scRNA@assays[["SCT"]]@SCTModel.list[["counts"]]@cell.attributes)
scRNA@assays[["SCT"]]@SCTModel.list[["counts"]]@model
###CLUSTER####

scRNA <- FindNeighbors(scRNA, dims = 1:10)
scRNA <- FindClusters(scRNA, resolution = 0.5)

?RunTSNE()

PrintTSNEParams

??FIt-SNE