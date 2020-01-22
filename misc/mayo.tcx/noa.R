library(EnsDb.Hsapiens.v79)
tbl.orig <- read.table("corrected_data.tsv", sep="\t", as.is=TRUE, nrow=-1)
dim(tbl.orig)
ensg <- rownames(tbl.orig)
ensg <- sub("\\..*$", "", ensg)
rownames(tbl.orig) <- ensg
tbl.map <- select(EnsDb.Hsapiens.v79, key=ensg,  columns=c("SYMBOL"),  keytype="GENEID")
dim(tbl.orig)
dim(tbl.map)

   #--------------------------------------------------------------------------------
   # not all of the ensg ids have gene symbols.  remove those from the tbl.orig
   #--------------------------------------------------------------------------------

map.failures <- setdiff(rownames(tbl.orig), tbl.map$GENEID)
map.failure.indices <- match(map.failures, rownames(tbl.orig))
tbl <- tbl.orig[-map.failure.indices,]
dim(tbl)
dim(tbl.map)

   #--------------------------------------------------------------------------------
   # all the ensg row names now have geneSymbol names, but some of the lattter
   # appear more than once.  rownames must be unique.  for instance, in tbl.map
   #               GENEID     SYMBOL
   # 2005 ENSG00000115239 GPR75-ASB3
   # 2006 ENSG00000270898 GPR75-ASB3
   # rows in tbl.map correspond to rows in tbl, so use which(duplicated(tbl.map$SYMBOL))
   # to eliminate rows in tbl
   # better yet: first sort on median or variance of expression of both genes
   #--------------------------------------------------------------------------------

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
   # looks like i got it wrong
   #-------------------------------------------------------------------------

stopifnot(all(mtx.from.noa >= 0))
minValue <- min(mtx.from.noa[mtx.from.noa > 0])

if(minValue <= 0){
  mtx.from.noa <- mtx.from.noa + abs(minValue) + .Machine$double.eps
  }

boxplot(as.numeric(mtx.from.noa), main="mtx.from.noa")
cor(mtx.from.noa[g1,], mtx.from.noa[g2,],method="spearman")   # 0.76  0.74
quartz()

   #----------------------------------------------------------------------------
   # max says the right approach is:  t(scale(t(mtx.log)))
   #----------------------------------------------------------------------------

mtx.log <- log10(mtx.from.noa)
quartz()

plot(mtx.from.noa[g1,], mtx.from.noa[g2,])   # 0.76
plot(mtx.scaled.log[g1,], mtx.scaled.log[g2,])

plot(mtx.from.noa[g1,], mtx.log[g1,])   # 0.76
plot(mtx.from.noa[g2,], mtx.log[g2,])   # 0.76
boxplot(as.numeric(mtx.log), main="mtx.log")
length(which(is.na(mtx.log)))
cor(mtx.log[g1,], mtx.log[g2,], method="spearman" )

mtx.scaled.log <- t(scale(t(mtx.log)))
quartz()
boxplot(as.numeric(mtx.scaled.log), main="mtx.scaled.log")
length(which(is.na(mtx.scaled.log))) # 982,500
cor(mtx.scaled.log[g1,], mtx.scaled.log[g2,], method="spearman")   # NA


