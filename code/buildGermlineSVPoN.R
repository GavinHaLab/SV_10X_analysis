#' buildLongRangerGermlineSVPoN.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:	  June 20, 2019
#' description: Build a panel of normal SV breakpoint list with frequencies in the matched normal cohort. Also, outputs binned regions with SV frequencies from the matched normal cohort that can be used as a black list.


library(optparse)
option_list <- list(
	make_option(c("--SVABAdir"), type="character", help="Path to directory containing SVABA SV results."),
	make_option(c("--LRdir"), type="character", help="Path to directory containing SVs extracted from LongRanger results."),
	make_option(c("--svaba_funcs"), type = "character", help = "Path to file containing SVABA R functions to source."),
	make_option(c("--blackListBinWidth"), type="integer", default=1000, help="Length (bp) for individual bins in the output black list. [Default: %default]"),
	make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze; string [Default: %default"),
	make_option(c("--outputPoNFile"), type="character", help="Path to output SV Panel of Normals (PoN) file."),
	make_option(c("--outputBlackListFile"), type="character", help="Path to output SV blacklist regions file.")
	)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(stringr)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

source(opt$svaba_funcs)


LRdir <- opt$LRdir
SVABAdir <- opt$SVABAdir
blackListBinWidth <- opt$blackListBinWidth
outputPoNFile <- opt$outputPoNFile
outputBlackListRegionsFile <- opt$outputBlackListFile
outImage <- paste0(outputPoNFile,".RData")
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
#chrs <- c(1:22, "X")
chrs <- as.character(eval(parse(text = opt$chrs)))
seqlevelsStyle(chrs) <- genomeStyle
bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
	seqinfo <- Seqinfo(genome=genomeBuild)
} else {
	seqinfo <- seqinfo(get(bsg))
}
seqinfo <- seqinfo[chrs]

buffer <- 100
minQual <- 20
minSVABA.span <- 0

save.image(outImage)

normLRFiles <- list.files(LRdir, pattern=".LR.germline.sv.txt", recursive = TRUE, full.names = TRUE)
names(normLRFiles) <- gsub(".LR.germline.sv.txt", "", basename(normLRFiles))
normSVABAFiles <- list.files(SVABAdir, pattern=".svaba.germline.sv.vcf", recursive = TRUE, full.names = TRUE)
names(normSVABAFiles) <- gsub(".svaba.germline.sv.vcf", "", basename(normSVABAFiles))

sampleList <- intersect(names(normLRFiles), names(normSVABAFiles))

sv.LR <- data.table()
sv.SVABA <- data.table()
for (i in 1:length(sampleList)){
	id <- sampleList[i]
	message("Loading LongRanger and SVABA files for ", id)
	# load and combine LR SV results
	svSample <- fread(normLRFiles[sampleList[i]])
	sv.LR <- rbind(sv.LR, cbind(Sample=id, svSample))
	# load and combine SVABA SV results
	svSample <- loadVCFtoDataTableByChromosome(normSVABAFiles[sampleList[i]], 
			chr=chrs, genomeStyle=genomeStyle, applyFilter = FALSE)
	sv.SVABA <- rbind(sv.SVABA, cbind(Sample = id, svSample))
}
save.image(outImage)
sv.LR <- unique(sv.LR)
#sv.SVABA <- unique(sv.SVABA)

sv.LR[, Tool := "LONGRANGER"]
sv.SVABA[, Tool := "SVABA"]
sv <- rbind(sv.LR, sv.SVABA[SPAN >= minSVABA.span | SPAN == -1], fill=TRUE)
sv[Tool == "SVABA", SOURCE := "SVABA"]
cols <- c("Sample", "chromosome_1", "start_1", "chromosome_2", "start_2")
cols.all <- c("Sample", "chromosome_1", "start_1", "chromosome_2", "start_2", "mateID", "FILTER", "Tool", "SOURCE", "DR", "SR", "AD", "PS", "SPAN", "orient_1", "orient_2")
sv <- sv[, cols.all, with = FALSE]

sv.1 <- sv[, !names(sv) %in% c("chromosome_2", "start_2"), with = FALSE]
sv.2 <- sv[, !names(sv) %in% c("chromosome_1", "start_1"), with = FALSE]
setnames(sv.1, c("chromosome_1", "start_1"), c("chromosome", "start"))
setnames(sv.2, c("chromosome_2", "start_2"), c("chromosome", "start"))
sv.1[, end := start]
sv.2[, end := start]
sv.dt <- rbind(sv.1, sv.2)

#sv.gr <- as(sv.dt, "GRanges")
#dup.ind <- which(duplicated(sv.gr) | duplicated(sv.gr, fromLast = TRUE)) 


counts <- sv.dt[, .(
	count = .N, SOURCE = paste0(unique(SOURCE), collapse=","),
	orient_1 = paste0(unique(orient_1), collapse=","),
	orient_2 = paste0(unique(orient_2), collapse=","),
	SPAN = paste0(sort(unique(SPAN)), collapse=",")
	), by=c("chromosome", "start")]
counts$chromosome <- factor(counts$chromosome, levels = chrs)
counts <- counts[order(chromosome, start)]

## get region overlap based on buffer region surrounding breakpoint 
#sv.dt.buff <- copy(sv.dt)
#sv.dt.buff[, start := spply(start - buffer, max, 0)]
#sv.dt.buff[, end := end + buffer]
sv.gr <- as(sv.dt, "GRanges")

tile.gr <- tileGenome(seqinfo, tilewidth=blackListBinWidth, cut.last.tile.in.chrom = TRUE)
hits.sv.count <- countOverlaps(query = tile.gr, subject = sv.gr)
tile.gr$SVbkptCounts <- hits.sv.count
hits <- as.data.table(findOverlaps(query = tile.gr, subject = sv.gr))
hits.sample.count <- hits[, sv.dt[subjectHits, length(unique(Sample))], by=queryHits]
tile.gr$SVsampleCounts <- 0
tile.gr$SVsampleCounts[hits.sample.count$queryHits] <- hits.sample.count$V1
tile.dt <- as.data.table(tile.gr)
# ## get all indices that are seen > 1 time(s)
# dup.ind <- which(duplicated(sv.gr) | duplicated(sv.gr, fromLast = TRUE)) 
save.image(outImage)

fwrite(counts, file = outputPoNFile, col.names=T, row.names=F, quote=F, sep="\t")
fwrite(tile.dt, file = outputBlackListRegionsFile, col.names=T, row.names=T, quote=F, sep="\t")



