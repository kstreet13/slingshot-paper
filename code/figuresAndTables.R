## ---- overall_setup ----
require(slingshot)
source('./code/helper_functions.R')


## ---- monocle_data_setup ----
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(HSMMSingleCell)
library(monocle)
data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))
marker_genes <- row.names(subset(fData(HSMM), 
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5", "ANPEP", "PDGFRA",
                                                        "MYOG", "TPM1", "TPM2", "MYH2", "MYH3", "NCAM1", "TNNT1", "TNNT2", "TNNC1",
                                                        "CDK1", "CDK2", "CCNB1", "CCNB2", "CCND1", "CCNA1", "ID1")))
diff_test_res <- differentialGeneTest(HSMM[marker_genes,], 
                                      fullModelFormulaStr="expression~Media")
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes <- merge(fData(HSMM), sig_genes, by="row.names")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.1))
ordering_genes <- intersect(ordering_genes, expressed_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, use_irlba=FALSE)
HSMM <- orderCells(HSMM, num_paths=2, reverse=TRUE)
HSMM <- HSMM[expressed_genes, pData(HSMM)$State != 3]
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(HSMM, fullModelFormulaStr="expression~Media")
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes <- merge(fData(HSMM), sig_genes, by="row.names")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.1))
ordering_genes <- intersect(ordering_genes, expressed_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, use_irlba=FALSE)
HSMM <- orderCells(HSMM, num_paths=2, reverse=TRUE)



## ---- waterfall_data_setup ----
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


## ---- multi-way_on_Monocle_data ----
X <- t(HSMM@reducedDimS); colnames(X) <- c('IC-1','IC-2')
lin <- get_lineages(X, pData(HSMM)$State)
crv <- get_curves(X, pData(HSMM)$State, lin)
plot_curves(X,pData(HSMM)$State,crv)

pst.mc <- monocle_pst(t(HSMM@reducedDimS))

pst.sl <- slingshot_pst(X,pData(HSMM)$State)

pairs(cbind(pData(HSMM)$Pseudotime,pst.mc,pst.sl),col=pData(HSMM)$State,pch=16)

pst.ts <- tscan_pst(all[,wfL1.idx], preproc = T)
df <- data.frame(Slingshot = pst.sl[wfL1.idx,1], Waterfall = pst.wf1, Monocle = pst.mc, TSCAN = pst.ts)
pairs(df, col=all.col[wfL1.idx], pch=16)




## ---- waterfall_1_lineage ----
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


## ---- waterfall_2_lineages


# ---- compare their two analyses to my one ----
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
crv <- get_curves(X, anno$V4, lin, extend = 'y',stretch = 2, shrink = F)
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



## ---- basic


## ---- Waterfall_data_PCA_figure ----
X <- prcomp(t(all))$x[,1:2]
lin <- get_lineages(X, anno$V4, start.clus = '1', end.clus = 'A')
crv <- get_curves(X, anno$V4, lin)
plot(X, col=all.col, pch=16, asp=1)
for(c in crv){lines(c$s[order(c$lambda),],col=1,lwd=2)}


## ---- multi-way_on_Waterfall_L1 ----
X <- prcomp(t(all))$x[,1:2]
pst.wf1 <- waterfall_pst(X[wfL1.idx,], k=4)
pst.sl <- slingshot_pst(X, anno$V4)
pst.mc <- monocle_pst(X[wfL1.idx,],num_paths = 1)
pst.ts <- tscan_pst(t(X[wfL1.idx,]))
df <- data.frame(Slingshot = pst.sl[wfL1.idx,1], Waterfall = pst.wf1, Monocle = pst.mc, TSCAN = pst.ts)
pairs(df, col=all.col[wfL1.idx], pch=16)

similarity <- apply(df,2,function(x){
  apply(df,2,function(y){
    Spp(x,y)
  })
})

# parallel lines
scaleAB <- function(x,a=0,b=1){
  ((x-min(x))/(max(x)-min(x)))*(b-a)+a
}
col.stabilize <- function(x){
  x <- 2*x
  y <- 27*x*x*(1-x)/4
  y[x > 2/3] <- 1
  return(y)
}
df01 <- data.frame(apply(df,2,scaleAB))
plot.new(); plot.window(xlim=0:1,ylim=c(.5,3.5))
abline(h=1:3,xlab='Pseudotime')
axis(1,lty=0); axis(2, at=1:3, labels = c('Waterfall','Slingshot','Monocle'),lty=0)
wdt <- matrix(scaleAB(c(abs(df01$Waterfall - df01$Slingshot), abs(df01$Monocle - df01$Slingshot)),1,3),ncol=2)
drk <- matrix(col.stabilize(c(abs(df01$Waterfall - df01$Slingshot), abs(df01$Monocle - df01$Slingshot))),ncol=2)
#drk <- matrix(scaleAB(c(abs(df01$Waterfall - df01$Slingshot), abs(df01$Monocle - df01$Slingshot))),ncol=2)
segments(df01$Waterfall,1, df01$Slingshot,2, col = alpha(1,drk[,1]), lwd=wdt[,1])
segments(df01$Monocle,3, df01$Slingshot,2, col = alpha(1,drk[,2]), lwd=wdt[,2])
points(df01$Waterfall,rep(1,sum(wfL1.idx)), pch=16, col=all.col[wfL1.idx], cex=1.5)
points(df01$Slingshot,rep(2,sum(wfL1.idx)), pch=16, col=all.col[wfL1.idx], cex=1.5)
points(df01$Monocle,rep(3,sum(wfL1.idx)), pch=16, col=all.col[wfL1.idx], cex=1.5)


## ---- 2_datasets_pairs ----
X1 <- t(HSMM@reducedDimS); colnames(X1) <- c('IC-1','IC-2')
clus1 <- pData(HSMM)$State
X2 <- prcomp(t(all))$x[wfL1.idx,1:2]
clus2 <- anno$V4[wfL1.idx]
  
pst.mc <- monocle_pst(X1)
pst.sl <- slingshot_pst(X1,clus1)
pst.wf <- waterfall_pst(X1,k=3)
pst.ts <- tscan_pst(exprs(HSMM), preproc = T)
df1 <- data.frame(Slingshot = pst.sl, Waterfall = pst.wf1, Monocle = pst.mc, TSCAN = pst.ts)

pst.mc <- monocle_pst(X2)
pst.sl <- slingshot_pst(X2,clus2)
pst.wf <- waterfall_pst(X2,k=4)
pst.ts <- tscan_pst(all[,wfL1.idx], preproc = T)
df2 <- data.frame(Slingshot = pst.sl[,1], Waterfall = pst.wf1, Monocle = pst.mc, TSCAN = pst.ts)


