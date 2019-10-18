#' annotateSVinPoN.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:	  July 15, 2019
#' description: Annotate SV events that overlap the Panel of Normals with SV frequencies observed in matched normal samples and black list regions. Filter SVs from input bedpe file. Filter based on cohort frequencies of SVs in PoN and SVs in bins (blacklist regions). An event (breakpoint pair) is filtered out of the results if 1 or both of the breakpoints is observed in the PoN SV breakpoint set or BlackList regions. PoN SV set and Blacklist regions to consider for filtering are determined based on the cohort frequency and whether it exceeds user-specified thresholds.


library(optparse)
option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--svaba_funcs"), type = "character", help = "Path to file containing SVABA R functions to source."),
	make_option(c("--svFile"), type="character", help = "Combined somatic SV calls for tumor sample."),
	make_option(c("--PoNFile"), type="character", help="Path to SV Panel of Normals (PoN) file."),
	make_option(c("--blackListFile"), type="character", help="Path to SV blacklist regions file."),
	make_option(c("--minFreqPoNSVBkptOverlap"), type="numeric", default=0.05, help="Minimum frequency of overlapping breakpoint across PoN cohort to consider for filtering out of results. Default [%default]"),
	make_option(c("--minFreqPoNBlackList"), type="numeric", default=0.75, help="Minimum frequency across PoN cohort to include blacklist region for filtering out of results. Default [%default]"),
	make_option(c("--outputSVAnnotFile"), type="character", help="Path to output SV file with PoN annotations."),
	make_option(c("--outputSVFiltFile"), type="character", help="Path to output SV file after filtering."),
	make_option(c("--outputSummary"), type="character", help="Path to output summary file.")
	)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

source(opt$svaba_funcs)

library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(stringr)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

id <- opt$id
svFile <- opt$svFile
ponFile <- opt$PoNFile
blackListFile <- opt$blackListFile


genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
minFreqPoNSVBkptOverlap <- opt$minFreqPoNSVBkptOverlap
minFreqPoNBlackList <- opt$minFreqPoNBlackList
outputSVAnnotFile <- opt$outputSVAnnotFile
outputSVFiltFile <- opt$outputSVFiltFile
outSummary <- opt$outputSummary

outImage <- paste0(outputSVFiltFile,".RData")
save.image(outImage)

############################################################
######### ANNOTATE PON SV AND BLACKLIST COUNTS #############
############################################################
sv <- fread(svFile)
pon <- fread(ponFile)
blist <- fread(blackListFile)

sv.1 <- sv[, c("chromosome_1", "start_1"), with = FALSE]
sv.2 <- sv[, c("chromosome_2", "start_2"), with = FALSE]
setnames(sv.1, c("chromosome_1", "start_1"), c("chromosome", "start"))
setnames(sv.2, c("chromosome_2", "start_2"), c("chromosome", "start"))
sv.1[, end := start]
sv.2[, end := start]
#sv.dt <- rbind(sv.1, sv.2)
sv.1.gr <- as(sv.1, "GRanges")
sv.2.gr <- as(sv.2, "GRanges")

## find SV breakpoints that are in the PoN ##
pon[, end := start]
hits1 <- findOverlaps(query = sv.1.gr, subject = as(pon, "GRanges"))
hits2 <- findOverlaps(query = sv.2.gr, subject = as(pon, "GRanges"))
sv.1[, SV.PoN.count := as.integer(0)]
sv.1[queryHits(hits1), SV.PoN.count := pon[subjectHits(hits1), count]]
sv.2[, SV.PoN.count := as.integer(0)]
sv.2[queryHits(hits2), SV.PoN.count := pon[subjectHits(hits2), count]]

## find number of SVs in blacklist regions ##
hits1 <- findOverlaps(query = sv.1.gr, as(blist, "GRanges"))
hits2 <- findOverlaps(query = sv.2.gr, as(blist, "GRanges"))
sv.1[, SV.blacklist.count := as.integer(0)]
sv.1[queryHits(hits1), SV.blacklist.count := blist[subjectHits(hits1), SVsampleCounts]]
sv.2[, SV.blacklist.count := as.integer(0)]
sv.2[queryHits(hits2), SV.blacklist.count := blist[subjectHits(hits2), SVsampleCounts]]
sv.1[, SV.blacklist.medianSVs := as.integer(0)]
sv.1[queryHits(hits1), SV.blacklist.count := median(blist[subjectHits(hits1), SVbkptCounts], na.rm=T)]
sv.2[, SV.blacklist.medianSVs := as.integer(0)]
sv.2[queryHits(hits2), SV.blacklist.count := median(blist[subjectHits(hits2), SVbkptCounts], na.rm=T)]
#tile.dt.bl <- merge(x=tile.dt, y=blist, by=c("seqnames", "start", "end", "width", "strand"), suffixes = c("", ".blacklist"))

sv[, SV.PoN.count_1 := sv.1$SV.PoN.count]
sv[, SV.PoN.count_2 := sv.2$SV.PoN.count]
sv[, SV.blacklist.count_1 := sv.1$SV.blacklist.count]
sv[, SV.blacklist.count_2 := sv.2$SV.blacklist.count]
sv[, SV.blacklist.medianSVs_1 := sv.1$SV.blacklist.medianSVs]
sv[, SV.blacklist.medianSVs_2 := sv.2$SV.blacklist.medianSVs]

#fwrite(sv, file = outputSVFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
writeBedpeToFile(sv, file=outputSVAnnotFile)

###########################################################
######## FILTER SV BASED ON PON AND BLACK LIST ############
###########################################################
numSamples <- length(unique(sv$Sample))
minNumPoNSV <- ceiling(numSamples * minFreqPoNSVBkptOverlap)
minNumPoNBlackList <- ceiling(numSamples * minFreqPoNBlackList)

sv.filt <- sv[!(start_1 == 0 | start_2 == 0)]

germSV.ind <- sv.filt[SV.PoN.count_1 >= minFreqPoNSVBkptOverlap | SV.PoN.count_2 >= minFreqPoNSVBkptOverlap | 
	SV.blacklist.count_1 >= minFreqPoNBlackList | SV.blacklist.count_2 >= minFreqPoNBlackList,
	which = TRUE]

## output filtered SV table to file ##
fwrite(sv.filt[-germSV.ind, ], file=outputSVFiltFile, col.names=T, row.names=F, quote=F, sep="\t")


## collect summary counts ##
germSV.tool.counts <- sv.filt[germSV.ind, table(Sample, Tool)]
germSV.support.counts <- sv.filt[germSV.ind, table(Sample, support)]
germSV.type.counts <- sv.filt[germSV.ind, table(Sample, CN_overlap_type)]
germSV.filter.counts <- sv.filt[germSV.ind, table(Sample, FILTER)]
summTab <- data.table(cbind(Sample=rownames(germSV.tool.counts), germSV.tool.counts, germSV.support.counts, germSV.type.counts, germSV.filter.counts))

## output summary counts to file ##
fwrite(summTab, file=outSummary, col.names=T, row.names=F, quote=F, sep="\t")

# save image
save.image(outImage)
