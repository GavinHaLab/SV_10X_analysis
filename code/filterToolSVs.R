
#' annotateSVinPoN.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' 
#' author: Anna Hoge <ahoge@fredhutch.org>
#''
#' date:	  Dec 26, 2020
#'
#' description: #effectively just removes only longranger calls for this dataset since minMoleculeLength=100kbps
#

library(optparse)
option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--svFile"), type="character", help = "Somatic SV calls for tumor sample with PoN germline SV counts annotated. Output of 'annotPoNSVs'."),
	make_option(c("--minFreqPoNSVBkptOverlap"), type="numeric", default=2, help="Minimum frequency of overlapping breakpoint in PoN cohort to consider for filtering out of results. Default [%default]"),
	# make_option(c("--minFreqPoNCNVBkptOverlap"), type="numeric", default=26, help="Minimum frequency of overlapping CNVs in PoN cohort to consider for filtering out of results. Default [%default]"),
	make_option(c("--minFreqPoNBlackList"), type="numeric", default=5, help="Minimum frequency across PoN cohort to include blacklist region for filtering out of results. Default [%default]"),
	make_option(c("--outputSVFile"), type="character", help="Path to output SV file with PoN annotations."),
	make_option(c("--outputSummary"), type="character", help="Path to output summary file.")
	)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(data.table)
setDTthreads(1)

id <- opt$id
svFile <- opt$svFile
minFreqPoNSVBkptOverlap <- opt$minFreqPoNSVBkptOverlap
# minFreqPoNCNVBkptOverlap <- opt$minFreqPoNCNVBkptOverlap
minFreqPoNBlackList <- opt$minFreqPoNBlackList
outFile <- opt$outputSVFile
outSummary <- opt$outputSummary

#parameter settings
minSPAN <- 10000
minNumTools <- 2
minMoleculeLength <- 100000
minLeftOverSV <- 10
# remove SVABA and LongRanger barcode-based SVs if molecule length < minMoleculeLength
removeBX <- FALSE
removeLR <- TRUE
# removeLRcnSize <- 250000


#begin run
#load SV file
sv <- fread(svFile, na.string = ".")
sv <- sv[, SV.id := NULL]
sv <- cbind(SV.id = 1:nrow(sv), sv)
sv[, Tool.multi := Tool] 
sv[!is.na(overlap.SVABA.id) & !is.na(overlap.GROCSVS.id) & is.na(overlap.LONGRANGER.id), Tool.multi := "SVABA,GROCSVS"]
sv[!is.na(overlap.SVABA.id) & is.na(overlap.GROCSVS.id) & !is.na(overlap.LONGRANGER.id), Tool.multi := "SVABA,LONGRANGER"]
sv[!is.na(overlap.SVABA.id) & !is.na(overlap.GROCSVS.id) & !is.na(overlap.LONGRANGER.id), Tool.multi := "SVABA,GROCSVS,LONGRANGER"]
sv[is.na(overlap.SVABA.id) & !is.na(overlap.GROCSVS.id) & !is.na(overlap.LONGRANGER.id), Tool.multi := "GROCSVS,LONGRANGER"]


#begin filtering

###########################################################
######## ANNOT SV BASED ON PON AND BLACK LIST ############
###########################################################
# numSamples <- length(unique(sv$Sample))
# minNumPoNSV <- ceiling(numSamples * minFreqPoNSVBkptOverlap)
# minNumPoNBlackList <- ceiling(numSamples * minFreqPoNBlackList)
sv <- sv[!(start_1 == 0 | start_2 == 0)]
sv[is.na(SV.PoN.sampleCount_1), SV.PoN.sampleCount_1 := 0]
sv[is.na(SV.PoN.sampleCount_2), SV.PoN.sampleCount_2 := 0]
germSV.ind <- sv[Tool %in% c("SVABA", "GROCSVS") &
				 (SV.PoN.sampleCount_1 >= minFreqPoNSVBkptOverlap | 
						SV.PoN.sampleCount_2 >= minFreqPoNSVBkptOverlap) | 
					  (SV.blacklist.sampleCount_1 >= minFreqPoNBlackList | 
						SV.blacklist.sampleCount_2 >= minFreqPoNBlackList), which = TRUE]
# germSV.ind <- sv[(Tool %in% c("SVABA", "GROCSVS") | 
# 					(Tool == "LONGRANGER" & !support %in% c("CN", "NoCNsupport"))) &
# 				 (SV.PoN.sampleCount_1 >= minFreqPoNSVBkptOverlap | 
# 						SV.PoN.sampleCount_2 >= minFreqPoNSVBkptOverlap) | 
# 					  (SV.blacklist.sampleCount_1 >= minFreqPoNBlackList | 
# 						SV.blacklist.sampleCount_2 >= minFreqPoNBlackList), which = TRUE]

# lr.germCNV.ind <- sv[Tool == "LONGRANGER" & support %in% c("CN", "NoCNsupport") &
# 					 ((SV.PoN.sampleCount_1 >= minFreqPoNCNVBkptOverlap | 
# 						SV.PoN.sampleCount_2 >= minFreqPoNCNVBkptOverlap) |
# 					 (SV.PoN.sampleCount_1 >= minFreqPoNSVBkptOverlap & 
# 						SV.PoN.sampleCount_2 >= minFreqPoNSVBkptOverlap)), which = TRUE]

## label PoN filter category ##
if (length(germSV.ind) > 0){
	sv[germSV.ind, SV.Filter := "PoN"]
	# sv[lr.germCNV.ind, SV.Filter := "PoN"]
}

###########################################################
######### ANNOT SV BASED ON TOOL OVERLAP #################
###########################################################
## exclude events not found in < minNumTools
samples <- unique(sv$Sample)
numSamples <- length(samples)
samplesFiltered <- data.table(Samples = samples, Tool.filtered = "")
toolExclude.ind <- NULL
for (i in 1:numSamples){
	
    sv.sample <- sv[Sample == samples[i]]
	numSV.sample <- nrow(sv.sample)
	
    # only apply to samples with Mean.Molecule.Length > minMoleculeLength
	molLen <- sv.sample[, unique(Mean.Molecule.Length)][1]
	toolInfo <- sv.sample[, .(overlap.SVABA.id, overlap.GROCSVS.id, overlap.LONGRANGER.id)]
	toolExclude.sample.ind <- which(rowSums(!is.na(toolInfo)) < minNumTools)
	numSVexclude <- length(toolExclude.sample.ind)
	
    # molecules larger than minMoleculeLength && number of SVS leftover is > minLeftOverSV
	if (molLen > minMoleculeLength && (numSV.sample - numSVexclude) >= minLeftOverSV){ 
		toolExclude.ind <- c(toolExclude.ind, sv.sample[toolExclude.sample.ind, SV.id])
		samplesFiltered[Samples == samples[i], Tool.filtered := "2Tools"]
	
    }else{ # molecules smaller than minMoleculeLength
		bx.ind <- NULL
		if (removeBX){
			bx.ind <- sv.sample[Tool == "SVABA" & grepl("BX", support), which = TRUE]
		}
		lr.ind <- NULL
		if (removeLR){
			# For LongRanger, only consider removing with support: 
			# For events that do no have overlap with CN segments (i.e. "NoCNsupport", "LRSV", "LOCAL_ASM")
			lr.ind <- sv.sample[Tool == "LONGRANGER" & !support %in% c("CN", "SV,CN"), 
				which = TRUE]			
		}
		# only remove events that are called by ONLY 1 tool
		removeTool.ind <- intersect(toolExclude.sample.ind, c(bx.ind, lr.ind))
		toolExclude.ind <- c(toolExclude.ind, sv.sample[removeTool.ind, SV.id])
		samplesFiltered[Samples == samples[i], Tool.filtered := "BX"]
	}
}
if (length(toolExclude.ind) > 0){
	sv[toolExclude.ind, SV.Filter := "ToolBX"]
}

###########################################################
#### CATEGORIZE SV BASED ON SIZE AND CN OVERLAP TYPE ######
###########################################################
sv[is.na(SV.Filter) & (is.na(CN_overlap_type) | CN_overlap_type == "Unknown-ShortSVwithCN"), 
	SV.Filter := "ShortSV"]
# exclude LONGRANGER calls with only CN support (i.e. not SVs) that don't overlap
# with other tools
sv[Tool.multi == "LONGRANGER" & support %in% c("CN", "NoCNsupport"), SV.Filter := "LRCN"]

###########################################################
#### KEEP SVs ######
###########################################################
# for SVs flagged as PoN or ShortSV but predicted multiple tools --> keep
#sv[Tool.multi == "SVABA,GROCSVS" | Tool.multi == "SVABA,LONGRANGER" | 
#	Tool.multi == "SVABA,GROCSVS,LONGRANGER" | Tool.multi == "GROCSVS,LONGRANGER", SV.Filter := "KEEP"]

# rest of the SVs with not Filter classification, then keep 
sv[is.na(SV.Filter), SV.Filter := "KEEP"]

## output filtered sv file 
fwrite(sv, file = outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", na = ".")

sv.counts <- cbind(sv[, .N, by=Sample],
	N.all.PoN.removed = sv[SV.Filter == "PoN", .N, by=Sample]$N, 
	N.all.Tool.removed = sv[SV.Filter == "ToolBX", .N, by=Sample]$N, 
	N.all.PoN.Tool.filtered = sv[SV.Filter != "PoN" & SV.Filter != "ToolBX", .N, by=Sample]$N, 
	N.10kb = sv[(SPAN == -1 | SPAN > minSPAN), .N, by=Sample]$N, 
	N.10kb.PoN.removed = sv[(SPAN == -1 | SPAN > minSPAN) & SV.Filter == "PoN", .N, by=Sample]$N, 
	N.10kb.Tool.removed = sv[(SPAN == -1 | SPAN > minSPAN) & SV.Filter == "ToolBX", .N, by=Sample]$N,
	N.10kb.ShortSV.removed = sv[(SPAN == -1 | SPAN > minSPAN) & SV.Filter == "ShortSV", .N, by=Sample]$N,
	N.10kb.PoN.Tool.filtered = sv[(SPAN == -1 | SPAN > minSPAN) & SV.Filter == "KEEP", .N, by=Sample]$N, 
	Molecule.length = sv[, unique(Mean.Molecule.Length)], 
	Tool.filtered.mode = samplesFiltered$Tool.filtered)

## print out summary
fwrite(sv.counts, file = outSummary, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


## split files into individual samples
# dir.create(outDir)
# for (i in 1:numSamples){
# 	id <- samples[i]
# 	sv.sample <- sv.filt[Sample == id]
# 	outBedFile <- paste0(outDir, "/", id, ".svaba-titan.PoNfilter.Toolfilter.bedpe")
# 	fwrite(sv.sample, file = outBedFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# }





