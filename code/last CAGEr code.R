library(BiocInstaller)
biocLite("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
library(CAGEr)


browseVignettes("CAGEr")
# analysis with wt2-4, ko1-4, aml1-3
# generate CAGEset (container) object
myCAGEset <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10", inputFiles = c("nanoCAGE-trimmed-02.bam", "nanoCAGE-trimmed-03.bam", "nanoCAGE-trimmed-04.bam", "nanoCAGE-trimmed-05.bam", "nanoCAGE-trimmed-06.bam", "nanoCAGE-trimmed-07.bam", "nanoCAGE-trimmed-08.bam", "nanoCAGE-trimmed-09.bam", "nanoCAGE-trimmed-10.bam", "nanoCAGE-trimmed-11.bam"), inputFilesType = "bam", sampleLabels = c("wt2", "wt3", "wt4", "ko1", "ko2", "ko3", "ko4", "aml1", "aml2", "aml3"))
myCAGEset
# reade in the data from bam files
getCTSS(myCAGEset)
# extract coordinates and counts of all TSS tags
ctss <- CTSStagCount(myCAGEset)
head(ctss)
sampleLabels(myCAGEset)
# check correlation between samples
corr.m <- plotCorrelation(myCAGEset, samples = "all", method = "pearson")
librarySizes(myCAGEset)
# merge samples
mergeSamples(myCAGEset, mergeIndex = c(1,1,1,2,2,2,2,3,3,3), mergedSampleLabels = c("wt","ko","aml"))
librarySizes(myCAGEset)
setColors(myCAGEset, colors=c("blue","red","forestgreen"))
sampleLabels(myCAGEset)
# check correlation after merging
corr.m <- plotCorrelation(myCAGEset, what = "CTSS", values = "raw", samples = "all", method = "pearson", tagCountThreshold = 2)
# extract coordinates and counts of all TSS tags after merging
ctss_merg <- CTSStagCount(myCAGEset)
head(ctss_merg)
# normalise counts (powerLaw, to 10Mil)
plotReverseCumulatives(myCAGEset, fitInRange = c(10, 20000), onePlot = TRUE)
normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(10, 20000), alpha = 1.05, T = 10*10^6)
# extract coordinates and normalised counts of all TSS tags
ctss_norm <- CTSSnormalizedTpm(myCAGEset)
head(ctss_norm)
corr.m <- plotCorrelation(myCAGEset, what = "CTSS", values = "normalized", samples = "all", method = "pearson", tagCountThreshold = 2)
# exporting normalised CAGE TSS reads to BedGraph file: 
exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = TRUE)
exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = FALSE)
# TSS clustering:
clusterCTSS(object=myCAGEset, threshold = 1.7, nrPassThreshold = 1, thresholdIsTpm = TRUE, method = "distclu", maxDist = 25, removeSingletons = TRUE, keepSingletonsAbove = 5, maxLength = 1000)
# cumul dist of TSS clusters & promoter width:
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
# do these again if needed, as clust parameters changed..
tc_WT <- tagClusters(myCAGEset, sample = "wt", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
tc_KO <- tagClusters(myCAGEset, sample = "ko", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
tc_AML <- tagClusters(myCAGEset, sample = "aml", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
write.table(tc_WT, file='tc_WT.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(tc_KO, file='tc_KO.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(tc_AML, file='tc_AML.tsv', row.names=FALSE, quote=FALSE, sep='\t')
# exporting interquantile widths:
exportToBed(object = myCAGEset, what = "tagClusters", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = TRUE, oneFile = TRUE)
exportToBed(object = myCAGEset, what = "tagClusters", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = TRUE, oneFile = FALSE)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", tpmThreshold = 1, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", tpmThreshold = 5, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", tpmThreshold = 10, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", tpmThreshold = 50, qLow = 0.1, qUp = 0.9)

# concensus promoters with tpmTreshold=5:
aggregateTagClusters(myCAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
consensusCl_5 <- consensusClusters(myCAGEset)
head(consensusCl_5)
write.table(consensusCl_5, file='consensusCl_5.tsv', row.names=FALSE, quote=FALSE, sep='\t')
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")
quantilePositions(myCAGEset, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9)
exportToBed(object = myCAGEset, what = "consensusClusters", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = FALSE, oneFile = TRUE)
exportToBed(object = myCAGEset, what = "consensusClusters", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = FALSE, oneFile = FALSE)
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", tpmThreshold = 5, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", tpmThreshold = 10, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", tpmThreshold = 50, qLow = 0.1, qUp = 0.9)
# concensus promoters with tpmTreshold=10:
aggregateTagClusters(myCAGEset, tpmThreshold = 10, qLow = 0.1, qUp = 0.9, maxDist = 100)
consensusCl_10 <- consensusClusters(myCAGEset)
head(consensusCl_10)
write.table(consensusCl_10, file='consensusCl_10.tsv', row.names=FALSE, quote=FALSE, sep='\t')
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")
quantilePositions(myCAGEset, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9)
exportToBed(object = myCAGEset, what = "consensusClusters", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = FALSE, oneFile = TRUE)
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", tpmThreshold = 5, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", tpmThreshold = 10, qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", tpmThreshold = 50, qLow = 0.1, qUp = 0.9)
# concensus promoters with tpmTreshold=50:
aggregateTagClusters(myCAGEset, tpmThreshold = 50, qLow = 0.1, qUp = 0.9, maxDist = 100)
#.. not done yet

# expression profiling wt-ko:
# across concensus clusters
getExpressionProfiles(myCAGEset, what = "consensusClusters", tpmThreshold = 10, nrPassThreshold = 1, method = "som", xDim = 3, yDim = 3)
plotExpressionProfiles(myCAGEset, what = "consensusClusters")
getExpressionProfiles(myCAGEset, what = "consensusClusters", tpmThreshold = 10, nrPassThreshold = 1, method = "kmeans", xDim = 9)
plotExpressionProfiles(myCAGEset, what = "consensusClusters")
KmeansGroups <- extractExpressionClass(myCAGEset, what = "consensusClusters", which = "all")
head(KmeansGroups)
write.table(KmeansGroups, file='KmeansGroups.tsv', row.names=FALSE, quote=FALSE, sep='\t')
exportToBed(object = myCAGEset, what = "consensusClusters", colorByExpressionProfile = TRUE, oneFile = TRUE)
# across TSS clusters
getExpressionProfiles(myCAGEset, what = "CTSS", tpmThreshold = 5, nrPassThreshold = 1, method = "som", xDim = 3, yDim = 3)
plotExpressionProfiles(myCAGEset, what = "CTSS")
getExpressionProfiles(myCAGEset, what = "CTSS", tpmThreshold = 5, nrPassThreshold = 1, method = "kmeans", xDim = 9)
plotExpressionProfiles(myCAGEset, what = "CTSS")
KmeansGroups <- extractExpressionClass(myCAGEset, what = "CTSS", which = "all")
head(KmeansGroups)
write.table(KmeansGroups, file='KmeansGroups.tsv', row.names=FALSE, quote=FALSE, sep='\t')
exportToBed(object = myCAGEset, what = "CTSS", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = TRUE, oneFile = TRUE)
# shifting promoters:
# wt vs. ko (norm)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters") # already done
scoreShift(myCAGEset, groupX = "wt", groupY = "ko", testKS = TRUE, useTpmKS = TRUE)
shiftPromoters_wtko_Norm <- getShiftingPromoters(myCAGEset, tpmThreshold = 0, scoreThreshold = 0, fdrThreshold = 0.5)
write.table(shiftPromoters_wtko_Norm, file='shiftPromoters_wtko_Norm.tsv', row.names=FALSE, quote=FALSE, sep='\t')
# wt vs. ko (raw)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters") # already done
scoreShift(myCAGEset, groupX = "wt", groupY = "ko", testKS = TRUE, useTpmKS = FALSE)
shiftPromoters_wtko_Raw <- getShiftingPromoters(myCAGEset, tpmThreshold = 0, scoreThreshold = 0, fdrThreshold = 0.5)
write.table(shiftPromoters_wtko_Raw, file='shiftPromoters_wtko_Raw.tsv', row.names=FALSE, quote=FALSE, sep='\t')
# ko vs. aml (norm)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters") # already done
scoreShift(myCAGEset, groupX = "ko", groupY = "aml", testKS = TRUE, useTpmKS = TRUE)
shiftPromoters_koaml_Norm <- getShiftingPromoters(myCAGEset, tpmThreshold = 0, scoreThreshold = 0, fdrThreshold = 0.5)
write.table(shiftPromoters_koaml_Norm, file='shiftPromoters_koaml_Norm.tsv', row.names=FALSE, quote=FALSE, sep='\t')
# ko vs. aml (raw)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters") # already done
scoreShift(myCAGEset, groupX = "ko", groupY = "aml", testKS = TRUE, useTpmKS = FALSE)
shiftPromoters_koaml_Raw <- getShiftingPromoters(myCAGEset, tpmThreshold = 0, scoreThreshold = 0, fdrThreshold = 0.5)
write.table(shiftPromoters_koaml_Raw, file='shiftPromoters_koaml_Raw.tsv', row.names=FALSE, quote=FALSE, sep='\t')

exportToBed(object = myCAGEset, what = "CTSS", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = TRUE, oneFile = TRUE)
exportToBed(object = myCAGEset, what = "CTSS", qLow = 0.1, qUp = 0.9, colorByExpressionProfile = TRUE, oneFile = FALSE)





