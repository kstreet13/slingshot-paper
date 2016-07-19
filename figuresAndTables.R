## ---- monocle_data_setup
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


## ---- waterfall_data_setup
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
col_sel <-c("#4882C3","#F26A6A","#13751B","#FF6A00","#E2CF00","#980B0B")
all.col <-rep("#BFBFBF30",length(branch))
all.col[which(anno$V4=="1")] <-"#4882C3" #1 blue
all.col[which(anno$V4=="2")] <-"#F26A6A" #2 salmon
all.col[which(anno$V4=="3")] <-"#13751B" #3 darkgreen
all.col[which(anno$V4=="4")] <-"#FF6A00" #4 orange
all.col[which(anno$V4=="5")] <-"#E2CF00" #5 yellow
all.col[which(anno$V4=="A")] <-"#980B0B" #Ap darkred
names(all.col) <-anno$V1
all.reduc <-prcomp(t(all))
redX <- all.reduc$x[,1:2]

## ---- ???
library(slingshot)
lin2 <- get_lineages(redX,anno$V4,start.clus = '1',end.clus='A')
crv2 <- get_curves(redX,anno$V4,lin2)

