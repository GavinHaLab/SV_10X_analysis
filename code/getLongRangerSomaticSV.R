#' combineSVABAandTITAN.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:	  October 2, 2018
#' description: Compare tumor and normal Long Ranger SVs to identify somatic events. Combines 


library(optparse)
option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--tenX_funcs"), type = "character", help = "Path to file containing 10X R functions to source."),
	make_option(c("--tumLargeSVFile"), type="character", help = "Long Ranger large SV calls for tumor sample (large_sv_calls.bedpe)"),
	make_option(c("--normLargeSVFile"), type="character", help = "Long Ranger large SV calls for normal sample (large_sv_calls.bedpe)"),
	make_option(c("--tumDeletionFile"), type="character", help = "Long Ranger deletion calls for tumor sample (dels.vcf.gz)"),
	make_option(c("--normDeletionFile"), type="character", help = "Long Ranger deletion calls for normal sample (dels.vcf.gz)"),
	make_option(c("--includeShortDeletions"), type="logical", default = TRUE, help = "true if Long Ranger deletion calls are to be included in output. [Default: %default]"),
	make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze; string [Default: %default"),
	make_option(c("--outDir"), type="character", help="Path to output directory."),
	make_option(c("--outputSVFile"), type="character", help="Path to output SV file with new annotations."),
	make_option(c("--outputNormSVFile"), type="character", help="Path to output matched normal SVs with new annotations.")
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

source(opt$tenX_funcs)

id <- opt$id
tumLargeSVFile <- opt$tumLargeSVFile
normLargeSVFile <- opt$normLargeSVFile
tumDeletionFile <- opt$tumDeletionFile
normDeletionFile <- opt$normDeletionFile
includeShortDeletions <- opt$includeShortDeletions
outputSVFile <- opt$outputSVFile
outputNormSVFile <- opt$outputNormSVFile
outDir <- opt$outDir
dir.create(outDir, recursive = TRUE)
outImage <- paste0(outDir, "/", id, ".RData")
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- as.character(eval(parse(text = opt$chrs)))
seqlevelsStyle(chrs) <- genomeStyle
#seqinfo <- Seqinfo(genome=genomeBuild)

buffer <- 1000
minDelLength <- 1000
minQual <- 20

save.image(outImage)

tum.sv <- getSVfromBEDPE(tumLargeSVFile, skip=1, genomeStyle = genomeStyle)
norm.sv <- getSVfromBEDPE(normLargeSVFile, skip=1, genomeStyle = genomeStyle)
if (includeShortDeletions){
	del.tum <- readVcf(tumDeletionFile, genome=genomeBuild)
	del.tum <- getSVfromCollapsedVCF.LR(del.tum, chrs=chrs, genomeStyle = genomeStyle)
	del.norm <- readVcf(normDeletionFile, genome=genomeBuild)
	del.norm <- getSVfromCollapsedVCF.LR(del.norm, chrs=chrs, genomeStyle = genomeStyle)
}

save.image(outImage)

## combine SVs with short Deletions if specified to do so
if (includeShortDeletions){
	tum.sv.del <- rbind(tum.sv, del.tum[FILTER=="PASS" & SPAN >= minDelLength & QUAL >= minQual], fill=TRUE)
	tum.sv.del <- tum.sv.del[, .(chromosome_1, start_1, chromosome_2, start_2, mateID, FILTER, SOURCE,
		HAP_ALLELIC_FRAC, ALLELIC_FRAC, DR, SR, PS, SPAN, orient_1, orient_2)]
	norm.sv.del <- rbind(norm.sv, del.norm[FILTER=="PASS" & SPAN >= minDelLength & QUAL >= minQual], fill=TRUE)
	norm.sv.del <- norm.sv.del[, .(chromosome_1, start_1, chromosome_2, start_2, mateID, FILTER, SOURCE,
		HAP_ALLELIC_FRAC, ALLELIC_FRAC, DR, SR, PS, SPAN, orient_1, orient_2)]
}else{
	tum.sv.del <- tum.sv[, .(chromosome_1, start_1, chromosome_2, start_2, mateID, 
		HAP_ALLELIC_FRAC, ALLELIC_FRAC, DR, SR, PS, SPAN, orient_1, orient_2)]
	tum.sv.del[, FILTER := "."]
	norm.sv.del <- norm.sv[, .(chromosome_1, start_1, chromosome_2, start_2, mateID, 
		HAP_ALLELIC_FRAC, ALLELIC_FRAC, DR, SR, PS, SPAN, orient_1, orient_2)]
	norm.sv.del[, FILTER := "."]
}

overlap <- getOverlapSV(tum.sv.del, norm.sv.del, buffer.x = buffer, buffer.y = buffer)

tum.sv.del <- cbind(Sample=id, SV.id = 1:nrow(tum.sv.del), tum.sv.del)
tum.sv.del[, SOMATIC := FALSE]
tum.sv.del[!SV.id %in% overlap$overlap.ind, SOMATIC := TRUE]

## sort by chromosome_1 and start_1
tum.sv.del <- tum.sv.del[!is.na(chromosome_1) & !is.na(chromosome_2), ]
tum.sv.del <- tum.sv.del[chromosome_1 %in% chrs & chromosome_2 %in% chrs, ]
tum.sv.del <- tum.sv.del[order(chromosome_1, start_1)]
tum.sv.del[, SV.id := 1:nrow(tum.sv.del)] # reassign SV.id

## exclude events without orientation ##
tum.sv.del <- tum.sv.del[!is.na(orient_1) & !is.na(orient_2)]

write.table(tum.sv.del, file = outputSVFile, col.names=T, row.names=F, quote=F, sep="\t")
write.table(norm.sv.del, file = outputNormSVFile, col.names=T, row.names=F, quote=F, sep="\t")
save.image(outImage)



