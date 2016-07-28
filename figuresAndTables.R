## ---- monocle_data_setup
{
library(HSMMSingleCell)
library(monocle)
data(HSMM_expr_matrix)
data(HSMM_sample_sheet)
data(HSMM_gene_annotation)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as(as.matrix(HSMM_expr_matrix), "sparseMatrix"),
                       phenoData = pd, featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
valid_cells <- row.names(subset(pData(HSMM), Mapped.Fragments > 1000000))
HSMM <- HSMM[,valid_cells]
}

## ---- waterfall_data_setup
{
all <- read.table('./data/Waterfall_TPM.txt')
anno <-read.table(file="./data/group0614WholeGenome4.csv",sep=",")
all <-all[,match(anno$V1,colnames(all))]
c <- cor(all, method="pearson"); d <- dist(t(all)); hr <- hclust(d, method = "ward", members=NULL)
branch=cutree(hr,k=6)
anno$V4 <-"0"
anno$V4[which(branch=="1")] <-"1"
anno$V4[which(branch=="2")] <-"2"
anno$V4[which(branch=="3")] <-"3"
anno$V4[which(branch=="4")] <-"4"
anno$V4[which(branch=="5")] <-"5"
anno$V4[which(branch=="6")] <-"A"
all.col <-rep("#BFBFBF30",length(branch))
all.col[which(anno$V4=="1")] <-"#4882C3" #1 blue
all.col[which(anno$V4=="2")] <-"#F26A6A" #2 salmon
all.col[which(anno$V4=="3")] <-"#13751B" #3 darkgreen
all.col[which(anno$V4=="4")] <-"#FF6A00" #4 orange
all.col[which(anno$V4=="5")] <-"#E2CF00" #5 yellow
all.col[which(anno$V4=="A")] <-"#980B0B" #Ap darkred
names(all.col) <-anno$V1
rm(c,d,hr,branch)
}

## ---- waterfall_1_lineage
{
source('~/Projects/slingshot_extras/WaterfallSupplement/Waterfall.R')
# interesting: their version of the MST plot
#mst.of.classification(all,k=6,col=paste0(all.col,""),seed=1)
anno.sub <- anno[-which(anno$V4=="A"),]
anno.sub <- anno.sub[anno.sub$V1 %notin% c("C11","C16","C25","C28","C13"),] # Removing S0; far left located; please see the Supplementary Methods II.1.(3) Additional trajectories
anno.sub <- anno.sub[anno.sub$V1 %notin% c("C48","C44","C138"),] # Removing outliers
all.sub <- all[,match(anno.sub$V1,colnames(all))]
all.col.sub <- all.col[match(anno.sub$V1,names(all.col))]
pca.sub <- prcomp(t(all.sub))
plot(pca.sub$x[,1:2],col=all.col.sub,pch=19,cex=2)
r <-kmeans(pca.sub$x[,1:2],centers=4) #they chose 4 subjectively
mst.of.classification(all.sub,k=4,col=all.col.sub,seed=1)
pseudotime.df <-pseudotimeprog.foo(all.sub,k=4,seed=1,color=all.col)
pseudotime.df <- pseudotime.df[match(colnames(all.sub),rownames(pseudotime.df)),]

lin <- get_lineages(pca.sub$x[,1:2], anno.sub$V4, start.clus = '1')
crv <- get_curves(pca.sub$x[,1:2], anno.sub$V4, lin)
plot_curves(pca.sub$x[,1:2], anno.sub$V4, crv, threeD = F)

plot(pseudotime.df$pseudotime, crv[[1]]$pseudotime, col=all.col.sub, pch=16)

pca <- prcomp(t(all))
lin.all <- get_lineages(pca$x[,1:2], anno$V4, start.clus = '1', end.clus = 'A')
crv.all <- get_curves(pca$x[,1:2], anno$V4, lin.all)
plot_curves(pca$x[,1:2], anno$V4, crv.all, threeD = F)
ps.overlap <- crv.all[[1]]$pseudotime[colnames(all) %in% colnames(all.sub)]

plot(pseudotime.df$pseudotime, ps.overlap, col=all.col.sub, pch=16)
}

## ---- 
{
library(slingshot)
lin <- get_lineages(pca.sub$x[,1:2],anno$V4,start.clus = '1',end.clus='A')
crv <- get_curves(redX,anno$V4,lin2)
}
