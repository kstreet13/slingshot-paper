## ---- overall_setup ----
{
  require(slingshot)
  source('./helper_functions.R')
}

## ---- monocle_data_setup
{
library(monocle)
data("HSMM_expr_matrix")
data("HSMM_sample_sheet")
data("HSMM_gene_annotation")
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as(as.matrix(HSMM_expr_matrix), "sparseMatrix"),
                       phenoData = pd, featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial.size())
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

valid_cells <- row.names(subset(pData(HSMM), Cells.in.Well == 1 & Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
HSMM <- HSMM[,valid_cells]

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(t(scale(t(L))))
# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")
}

## ---- waterfall_data_setup
{
all <- read.table('./data/Waterfall_TPM.txt')
anno <-read.table(file="./data/group0614WholeGenome4.csv",sep=",")
source('../slingshot_extras/WaterfallSupplement/Waterfall.R')
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
wfALL.idx <- ! anno$V1 %in% c("C11","C16","C25","C28","C13","C48","C44","C138")
all <- all[,wfALL.idx]
anno <- anno[wfALL.idx,]
all.col <- all.col[wfALL.idx]
wfL1.idx <- (anno$V4 != 'A')
wfL2.idx <- (anno$V4 %in% c('1','2','3','A'))
}

## ---- waterfall_1_lineage
{
# compare to published analysis of first lineage
X <- prcomp(t(all)[wfL1.idx,])$x[,1:2]
pst.wf <- waterfall_pst(X, k=4) #they chose 4 (subjectively), 5 looks pretty similar
pst.sl <- slingshot_pst(X, anno$V4[wfL1.idx])
plot(pst.wf, pst.sl, col=all.col[wfL1.idx], pch=16)
Spp(pst.wf,pst.sl)

# compare to analogous analysis of second lineage
X <- prcomp(t(all)[wfL2.idx,])$x[,1:2]
pst.wf <- waterfall_pst(X, k=4) # 4 makes sense
pst.sl <- slingshot_pst(X, anno$V4[wfL2.idx])
plot(pst.wf, pst.sl, col=all.col[wfL2.idx], pch=16)
Spp(pst.wf,pst.sl)
}

## ---- waterfall_2_lineages
{
# compare their two analyses to my one
X <- prcomp(t(all))$x[,1:2]
pst.wf1 <- waterfall_pst(X[wfL1.idx,], k=4)
pst.wf2 <- waterfall_pst(X[wfL2.idx,], k=4)
pst.sl <- slingshot_pst(X, anno$V4)
plot(pst.wf1, pst.sl[wfL1.idx,1], col=all.col[wfL1.idx], pch=16)
plot(pst.wf2, pst.sl[wfL2.idx,2], col=all.col[wfL2.idx], pch=16)
Spp(pst.wf1,pst.sl[wfL1.idx,1])
Spp(pst.wf2,pst.sl[wfL2.idx,2])

# pca plot with tree, smooth curves for figure?
lin <- get_lineages(X, anno$V4)
crv <- get_curves(X, anno$V4, lin, extend = 'y',stretch = 0)
forest <- lin$forest
nclus <- nrow(forest)
centers <- t(sapply(rownames(forest),function(clID){
  x.sub <- X[anno$V4 == clID,]
  return(colMeans(x.sub))
}))
center.col <- sapply(rownames(forest),function(clID){
  all.col[which.max(anno$V4 == clID)]
})
plot(X, col=all.col, cex=2, pch=16, asp = 1)
for(i in 1:(nclus-1)){
  for(j in (i+1):nclus){
    if(forest[i,j]==1){
      lines(centers[c(i,j),1], centers[c(i,j),2], lwd=3)
    }
  }
}
#points(centers, cex = 3, col = center.col, pch=16)
points(centers, cex = 2, pch=16)
#text(centers, labels = rownames(centers), adj = c(1,1))
for(l in 1:length(crv)){ lines(crv[[l]]$s, lwd=3, col=2)}
}

## ---- 
{
library(slingshot)
lin <- get_lineages(pca.sub$x[,1:2],anno$V4,start.clus = '1',end.clus='A')
crv <- get_curves(redX,anno$V4,lin2)
}
