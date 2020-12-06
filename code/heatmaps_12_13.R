--12--
# 1 heatmap for DE pattern identification
blah <- read.delim('cagen_log_rpkm.txt', header=TRUE, stringsAsFactors=FALSE)
all.de <- list(blah$Symbol[1:162])

every.g <- blah[,2:7]
rownames(every.g) <- blah$Symbol

heat.ab <- every.g - rowMeans(every.g)

# can try doing with either original (every.g) or centered (heat.ab)

rowseps <- cumsum(sapply(all.de, FUN=length))
last <- 0L
re.heat <- list()
for (x in all.de) {
    current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
    if (length(x)!=1L) {
      ddr <- as.dendrogram(hclust(dist(current), method="average"))
      ddr <- reorder(ddr, rowMeans(current))
      current <- current[order.dendrogram(ddr),,drop=FALSE] 
    }
    re.heat[[length(re.heat)+1L]] <- current
    last <- last + length(x) 
    }
heat.ab <- do.call(rbind, re.heat)

lwid <- c(4, 0.3*ncol(heat.ab) + 2)
lhei <- c(2, 0.2*nrow(heat.ab) + 2)

cur.fname <- "cagenes_heat_av_dendrpo.pdf"
pdf(cur.fname, height=sum(lhei), width=sum(lwid))

cur.fname <- "cagenes_heat.png"
png(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

cur.fname <- "cagenes_heatt.tiff"
tiff(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

require(gplots)
heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
heatmap.2(heat.abb, dendrogram="row", trace="none", density.info="none",
          col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
          rowsep=rowseps, Colv=FALSE, Rowv=T, sepwidth=c(0.2, 0.2), margin=c(10,10))
dev.off()

# 2 Now do the separate heatmaps for clusters 1, 2+6, 5 (DNclust); 3, 4+8, 7 (UPclust); 0 (0Clust):
blah <- read.delim('cagen_UP_DN_cl.0.txt', header=TRUE, stringsAsFactors=FALSE)
all.de <- list(blah$Symbol)

every.g <- blah[,2:7]
rownames(every.g) <- blah$Symbol

heat.ab <- every.g - rowMeans(every.g)

rowseps <- cumsum(sapply(all.de, FUN=length))
last <- 0L
re.heat <- list()
for (x in all.de) {
  current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
  if (length(x)!=1L) {
    ddr <- as.dendrogram(hclust(dist(current), method = "complete"))
    ddr <- reorder(ddr, rowMeans(current))
    # ddr <- reorder(ddr, rowMeans(current[,5:6,drop=FALSE])-rowMeans(current[,1:4,drop=FALSE]))
    current <- current[order.dendrogram(ddr),,drop=FALSE] 
  }
  re.heat[[length(re.heat)+1L]] <- current
  last <- last + length(x) 
}
heat.ab <- do.call(rbind, re.heat)

lwid <- c(4, 0.3*ncol(heat.ab) + 2)
lhei <- c(2, 0.2*nrow(heat.ab) + 2)

cur.fname <- "cagen_heat_cl0.pdf"
pdf(cur.fname, height=sum(lhei), width=sum(lwid))

cur.fname <- "cagen_heat_cl0.png"
png(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

cur.fname <- "cagen_heat_cl0.tiff"
tiff(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

require(gplots)
heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
          col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
          rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10))
dev.off()
rm(x)

# 3 Now do the separate heatmaps for Hallmarks:
blah <- read.delim('H47_log_rpkm.txt', header=TRUE, stringsAsFactors=FALSE)
all.de <- list(blah$Symbol) # can change it to $GeneName if required

every.g <- blah[,2:7]
rownames(every.g) <- blah$Symbol # can change it to $GeneName if required

heat.ab <- every.g - rowMeans(every.g)

rowseps <- cumsum(sapply(all.de, FUN=length))
last <- 0L
re.heat <- list()
for (x in all.de) {
  current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
  if (length(x)!=1L) {
    ddr <- as.dendrogram(hclust(dist(current), method = "complete"))
    ddr <- reorder(ddr, rowMeans(current[,1:2,drop=FALSE])-rowMeans(current))
    #ddr <- reorder(ddr, rowMeans(current[,5:6,drop=FALSE])-rowMeans(current[,1:4,drop=FALSE]))
    current <- current[order.dendrogram(ddr),,drop=FALSE] 
  }
  re.heat[[length(re.heat)+1L]] <- current
  last <- last + length(x) 
}
heat.ab <- do.call(rbind, re.heat)
heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)

lwid <- c(4, 0.3*ncol(heat.ab) + 2)
lhei <- c(2, 0.2*nrow(heat.ab) + 2)

cur.fname <- "H47_heat.pdf"
pdf(cur.fname, height=sum(lhei), width=sum(lwid))

cur.fname <- "H47_heat.png"
png(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

cur.fname <- "H47_heat.tiff"
tiff(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

require(gplots)
heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
          col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
          rowsep=rowseps, Colv=FALSE, Rowv=F, sepwidth=c(0.2, 0.2), margin=c(10,10))
dev.off()
rm(x)

###

set.seed(101)
kClust <- kmeans(heat.ab, 6, nstart = 50)
kClust$size
kClust$cluster
write.table(as.data.frame(kClust$cluster), file='H10_kClust.txt', row.names=T, sep='\t', quote=FALSE)



--13--
blah <- read.delim('H_pathInCa.txt', header=TRUE, stringsAsFactors=FALSE)
all.de <- list(blah$Symbol) # can change it to $GeneName if required

every.g <- blah[,2:7]
rownames(every.g) <- blah$Symbol # can change it to $GeneName if required

heat.ab <- every.g - rowMeans(every.g)

rowseps <- cumsum(sapply(all.de, FUN=length))
last <- 0L
re.heat <- list()
for (x in all.de) {
  current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
  if (length(x)!=1L) {
    # ddr <- as.dendogram(hcust(as.dist(1-cor(t(current))), method = "average")) # try this instead !
    ddr <- as.dendrogram(hclust(dist(current), method = "complete"))
    #ddr <- reorder(ddr, rowMeans(current))
    ddr <- reorder(ddr, -rowMeans(current[,1:2,drop=FALSE])-rowMeans(current))
    #ddr <- reorder(ddr, rowMeans(current[,1:2,drop=FALSE])-rowMeans(current[,4:6,drop=FALSE]))
    current <- current[order.dendrogram(ddr),,drop=FALSE] 
  }
  re.heat[[length(re.heat)+1L]] <- current
  last <- last + length(x) 
}
heat.ab <- do.call(rbind, re.heat)
heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)

lwid <- c(4, 0.3*ncol(heat.ab) + 2)
lhei <- c(2, 0.2*nrow(heat.ab) + 2)

cur.fname <- "HD_pathInCa.pdf"
pdf(cur.fname, height=sum(lhei), width=sum(lwid))

setEPS()
cur.fname <- "HD_pathInCa.pdf"
postscript(cur.fname, onefile = FALSE, width = sum(lwid), height = sum(lhei), paper = "special")
# with ggplot2 this would work: ggsave(file="name.eps") 

cur.fname <- "HD_pathInCa.png"
png(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

cur.fname <- "HD_pathInCa.tiff"
tiff(cur.fname, height=sum(lhei), width=sum(lwid), units="in", res=300)

require(gplots)
heatmap.2(heat.abb, dendrogram="row", trace="none", density.info="none",
          col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
          rowsep=rowseps, Colv=FALSE, Rowv=T, sepwidth=c(0.2, 0.2), margin=c(10,10))
dev.off()

#clusterCut <- cutree(current, 10)
#write.table(as.data.frame(clusterCut), file='PK_Clust.txt', row.names=T, sep='\t', quote=FALSE)

rm(x)

