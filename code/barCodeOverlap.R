#' barcodeOverlapFromVCF_v2.1.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:	  July 26, 2018
#' description: Performs barcode overlap of all pairs of breakpoints from 10X Genomics WGS data

library(optparse)
option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--tenX_funcs"), type = "character", help = "Path to file containing 10X R functions to source."),
	make_option(c("--svaba_funcs"), type = "character", help = "Path to file containing SVABA R functions to source."),
	make_option(c("--tumBam"), type="character", help = "Path to tumor LongRanger bam file"),
	make_option(c("--vcf"), type="character", help = "Path to unfiltered somatic sv vcf file. Should have suffix svaba.unfiltered.somatic.sv.vcf"),
	make_option(c("--bps"), type="character", help = "Path to full breakpoint file. Should have suffix bps.txt.gz"),
	make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze; string [Default: %default"),
	make_option(c("--minMapQ"), type="integer", default=20, help = "Minimum mapping quality to use for analysis. Default [%default]"),
	make_option(c("--minLength"), type="integer", default=10000, help = "Minimum length to consider for barcode rescue. Default [%default]"),
	make_option(c("--windowSize"), type="integer", default=1000, help = "Window size at each breakpoint end of the pair for an sv. Default [%default]"),
	make_option(c("--minReadOverlapSupport"), type="integer", default=2, help = "Minimum number of read overlap support to use at each breakpoint window for breakpoint pairs of an sv event. Default [%default]"),
	make_option(c("--outFile"), type="character", help = "Output vcf file.")
)
options(stringsAsFactors = FALSE, scipen = 999, width=175)

library(data.table)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(reshape2)
library(VariantAnnotation)
library(tools)
library(GenomeInfoDb)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

source(opt$svaba_funcs)
source(opt$tenX_funcs)

id <- opt$id
tumBamFile <- opt$tumBam
vcfFile <- opt$vcf
bpsFile <- opt$bps
minMAPQ <- opt$minMapQ
minLength <- opt$minLength
windowSize <- opt$windowSize
minReadOverlapSupport <- opt$minReadOverlapSupport
outFile <- opt$outFile
build <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- eval(parse(text = opt$chrs))
seqlevelsStyle(chrs) <- genomeStyle
outImage <- gsub("vcf", "RData", outFile)
save.image(outImage)

## load vcf file ##
vcf <- readVcf(vcfFile, genome = build)

## filter vcf by span length ##
indSPAN <- info(vcf)$SPAN >= minLength | info(vcf)$SPAN == -1 | fixed(vcf)$FILTER == "PASS"
## filter vcf by FILTER == DUPREADS and LOCALMATCH
indFILTER <- !rowRanges(vcf)$FILTER %in% c("DUPREADS", "LOCALMATCH")
## set vcf file to desired genomeStyle
vcf <- setVCFgenomeStyle(vcf, genomeStyle = genomeStyle)
## filter by autosomes and sex chromosomes
indCHR <- seqnames(rowRanges(vcf)) %in% chrs
ind <- indSPAN & indFILTER & indCHR
vcf <- vcf[ind]
message("Filtered ", sum(!ind), " rows - ", sum(ind), " remaining")
sv <- getSVfromCollapsedVCF(vcf, chrs = chrs, genomeStyle = genomeStyle)
## filter rows with NA in alt_1, alt_2, orient_1, orient_2 fields
sv <- sv[!is.na(alt_1) & !is.na(alt_2) & !is.na(orient_1) & !is.na(orient_2)]
setkey(sv, mateID)

## load bps.txt file ##
bps <- loadBPStoDataTableByChromosome(bpsFile, tumor.id = id, chrs = chrs, 
    minLength = minLength, dupSV.bpDiff = 1000)
bps[, mateID := 1:nrow(bps)]
bps[bxtags == "x", bxtags := NA]


## annotate bxtags from bps to sv
overlap <- annotVCFbyOverlap(sv, bps, annotCol="bxtags", buffer = 0)
sv[overlap$ind, bxtags := overlap$values]

## set up Rsamtools scan bam params
tags <- c("BX")
flags <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
        hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
        isFirstMateRead = NA, isSecondMateRead = NA, 
        isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
        isDuplicate = FALSE)
fields <- scanBamWhat()       

save.image(outImage)

## Get barcode overlap - main function call ## 
overlap.count <- getBXoverlap(sv, tumBamFile, minReadOverlapSupport = minReadOverlapSupport, windowSize = windowSize, flags = flags, fields = fields, tags = tags)
save.image(outImage)

## set up additional info field for new VCF ##
message("Adding BXOL to header")
vcf <- addInfoFieldtoCollapsedVCF(vcf, field.name="BXC.1", 
	field.desc="Number of barcodes from proper pairs at first breakpoint.", 
	field.defaultValue=0, field.type="Integer")
vcf <- addInfoFieldtoCollapsedVCF(vcf, field.name="BXC.2", 
	field.desc="Number of barcode from proper pairs at second breakpoint.", 
	field.defaultValue=0, field.type="Integer")
vcf <- addInfoFieldtoCollapsedVCF(vcf, field.name="BXOL", 
	field.desc="Number of barcode overlaps between proper pairs spanning both breakpoints.", 
	field.defaultValue=0, field.type="Integer")
vcf <- addInfoFieldtoCollapsedVCF(vcf, field.name="BXOL.plus.Support", 
	field.desc="Number of barcode overlaps between proper pairs, split reads, and discordant reads spanning both breakpoints.", 
	field.defaultValue=0, field.type="Integer")
vcf <- addInfoFieldtoCollapsedVCF(vcf, field.name="BXSupport", 
	field.desc="Number of barcode of SR and DR supporting both breakpoints.", 
	field.defaultValue=0, field.type="Integer")
vcf <- addInfoFieldtoCollapsedVCF(vcf, field.name="BXTags", 
	field.desc="Barcodes of SR and DR supporting both breakpoints.", 
	field.defaultValue=0, field.type="String", Number = ".")

## assign barcode overlap count to vcf ##
mateID1 <- gsub(":2", ":1", overlap.count$mateID) ## rownames only contains :2
info(vcf)[mateID1, "BXC.1"] <- overlap.count$Region1.BXcount
info(vcf)[mateID1, "BXC.2"] <- overlap.count$Region2.BXcount
info(vcf)[mateID1, "BXOL"] <- overlap.count$OverlapCount.bxol.only
info(vcf)[mateID1, "BXOL.plus.Support"] <- overlap.count$OverlapCount.bxol.plus.support
info(vcf)[mateID1, "BXSupport"] <- overlap.count$Support.BXCount
info(vcf)[mateID1, "BXTags"] <- sv[as.character(overlap.count$mateID), bxtags]
mateID2 <- overlap.count$mateID ## use :1
info(vcf)[mateID2, "BXC.1"] <- overlap.count$Region1.BXcount
info(vcf)[mateID2, "BXC.2"] <- overlap.count$Region2.BXcount
info(vcf)[mateID2, "BXOL"] <- overlap.count$OverlapCount.bxol.only
info(vcf)[mateID2, "BXOL.plus.Support"] <- overlap.count$OverlapCount.bxol.plus.support
info(vcf)[mateID2, "BXSupport"] <- overlap.count$Support.BXCount
info(vcf)[mateID2, "BXTags"] <- sv[as.character(overlap.count$mateID), bxtags]

sv <- getSVfromCollapsedVCF(vcf, genomeStyle = genomeStyle)
save.image(outImage)

## output to VCF file ##
writeVcf(vcf, filename = outFile)

## output sv object as a table ##
outFileSV <- gsub("vcf", "txt", outFile)
write.table(sv, file=outFileSV, col.names=T, row.names=F, quote=F, sep="\t")

#mat <- as.data.frame(cbind(starts = starts[-1], counts = overlap.count[, 3]))
#mat <- mat[!is.na(mat[, 2]), ]
#gp <- ggplot(mat, aes(x=starts, y=counts)) + geom_point() + geom_line() +
#						ylab("Barcode Overlap Count") + xlab("ChrX Coordinates") +
#						theme_bw()
#outPlot <- paste0(outPrefix, "_", eventSize * windowSize/1000, "kb.pdf")
#ggsave(outPlot, width = 10, height = 3, units="in")
	
	
	
