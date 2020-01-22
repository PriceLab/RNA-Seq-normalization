library(EnsDb.Hsapiens.v79)
tbl.orig <- read.table("corrected_data.tsv", sep="\t", as.is=TRUE, nrow=-1)
dim(tbl.orig)
ensg <- rownames(tbl.orig)
ensg <- sub("\\..*$", "", ensg)
rownames(tbl.orig) <- ensg
tbl.map <- select(EnsDb.Hsapiens.v79, key=ensg,  columns=c("SYMBOL"),  keytype="GENEID")
dim(tbl.orig)
dim(tbl.map)
map.failures <- setdiff(rownames(tbl.orig), tbl.map$GENEID)
map.failure.indices <- match(map.failures, rownames(tbl.orig))
tbl <- tbl.orig[-map.failure.indices,]
dim(tbl)
dups <- which(duplicated(tbl.map$SYMBOL))
length(dups)
tbl <- tbl[-dups,]
tbl.map <- tbl.map[-dups,]
dim(tbl); dim(tbl.map)
nrow(tbl.map)
indices <- match(rownames(tbl), tbl.map$GENEID)
length(indices)
dim(tbl)
rownames(tbl) <- tbl.map$SYMBOL[indices]
tbl[1:10, 1:10]
mtx.from.noa <- as.matrix(tbl)

g1 <- "MEF2C"
g2 <- "PKNOX2"

cor(mtx.from.noa[g1,], mtx.from.noa[g2,])   # 0.760881


   #------------------------------------------------------------------------
   # apply my understanding of max & michael's transformation of count data
   #-------------------------------------------------------------------------

stopifnot(all(mtx.from.noa >= 0))
minValue <- min(mtx.from.noa[mtx.from.noa > 0])

if(minValue <= 0){
  mtx.from.noa <- mtx.from.noa + abs(minValue) + .Machine$double.eps
  }

mtx.scaled <-  t(scale(t(mtx.from.noa)))
plot(mtx.scaled[g1,], mtx.scaled[g2,])
cor(mtx.scaled[g1,], mtx.scaled[g2,])   # 0.76

mtx.log.scaled <- log10(mtx.scaled)
plot(mtx.log.scaled[g1,], mtx.log.scaled[g2,])

mtx.log <- log10(mtx.from.noa)
plot(mtx.log[g1,], mtx.log[g2,])

mtx.scaled.log <- t(scale(t(mtx.log)))
plot(mtx.scaled.log[g1,], mtx.scaled.log[g2,])

mtx.scaled <-  t(scale(t(mtx.2)))
mtx.log.scaled <- t(scale(t(mtx.2)))

cor(mtx.log.scaled[g1,], mtx.log.scaled[g2,])   # 0.407936

mtx.old <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/temporalCortex.15167x264.RData"))
cor(mtx.old[g1,], mtx.old[g2,])   # 0.8315445

