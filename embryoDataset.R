counts <- read.table('~/Downloads/counts.txt')
gfilter <- apply(counts,1,function(x){
  sum(x >= 20) >= 20
})
counts <- counts[gfilter,]
meta <- read.table('~/Downloads/metadata.txt',header = TRUE)
require(RColorBrewer)
cc <- c(brewer.pal(9, "Set1")[-1], brewer.pal(8, "Set2"), brewer.pal(9, "Set3"), brewer.pal(6, "Dark2"))
batch <- sapply(as.character(meta[,1]),function(s){unlist(strsplit(s,split='[.]'))[2]})
meta$batch <- factor(batch); rm(batch)

good <- meta$cluster != 'grey'
depth <- log1p(colSums(counts))

require(scone)
#norm <- DESEQ_FN(counts)
norm <- FQ_FN(counts)

data("housekeeping")
gnames <- read.csv('~/Downloads/mart_export.txt',header = T)
gnames[,2] <- tolower(gnames[,2])
housekeeping <- tolower(housekeeping[,1])
hk <- as.character(gnames[gnames[,2] %in% housekeeping,1])
hk <- hk[hk %in% rownames(norm)]

require(RUVSeq)
ruv <- RUVg(round(as.matrix(norm)), hk, 1)

X <- prcomp(t(log1p(ruv$normalizedCounts)))$x
require(fastICA)
ica <- fastICA(t(log1p(ruv$normalizedCounts)),n.comp = 3)
X <- prcomp(t(log1p(norm)))$x

require(rgl)
plot3d(X[,1:3], col=as.character(meta$cluster), size = 3)


