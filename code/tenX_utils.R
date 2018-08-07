#' author: Gavin Ha 
#' 		Dana-Farber Cancer Institute
#'		  Broad Institute
#' contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
#' date:	  March 17, 2017

#' @import data.table
#' @import GenomicRanges

loadBXcountsFromBEDDir <- function(bxDir, chrs = c(1:22, "X", "Y"), minReads = 2){
	files <- list.files(bxDir, pattern=".bed", full.names = TRUE)
	message("Loading BX counts from ", bxDir)
	bxAll <- NULL
	for (i in files){
		message(i)
		awkcmd <- paste0("awk -F \"\t\" \'{if (!$4) {print $1\"\t\"$2\"\t\"$3\"\t\"\"NA\"} else {print $1\"\t\"$2\"\t\"$3\"\t\"$4}}\' ", i)
		bxChr <- fread(awkcmd, header = FALSE, na.strings = c("NA"))
		colnames(bxChr) <- c("chr", "start", "end", "BX")
		bxChr[, BXcounts:=sapply(BX, function(x){
		if (!is.na(x)){
				bxCodes <- unlist(tstrsplit(x,","))
				bxReads <- as.numeric(tstrsplit(bxCodes, "_")[[2]])
				return(sum(bxReads >= minReads))
			}else{
				return(0)
			}
		})]
		bxChr[, BX:=NULL]
		bxAll <- rbind(bxAll, bxChr)
	}
	bxGR <- RangedData(space = bxAll[[1]], ranges = IRanges(start = bxAll[[2]], 
										 end = bxAll[[3]]), BX.count = bxAll[[4]])
	bxGR <- keepChr(bxGR, chr = chrs)
	return(bxGR)
}


loadReadCountsFromBed <- function(counts, chrs = c(1:22, "X", "Y"), gc = NULL, map = NULL, centromere = NULL, flankLength = 100000, targetedSequences = NULL, genomeStyle = "NCBI"){
	require(HMMcopy)
	require(GenomeInfoDb)
	names(counts) <- setGenomeStyle(names(counts), genomeStyle)
	counts <- keepChr(counts, chrs)
	if (!is.null(gc)){ 
		names(gc) <- setGenomeStyle(names(gc), genomeStyle)
		gc <- keepChr(gc, chrs)
		counts$gc <- gc$value
	}
	if (!is.null(map)){ 
		names(map) <- setGenomeStyle(names(map), genomeStyle)
		map <- keepChr(map, chrs)
		counts$map <- map$value
	}
	colnames(counts)[1] <- c("reads")
	
	# remove centromeres
	if (!is.null(centromere)){ 
		counts <- excludeCentromere(counts, centromere, flankLength = flankLength)
	}
	# keep targeted sequences
	if (!is.null(targetedSequences)){
		countsExons <- filterByTargetedSequences(counts, targetedSequences)
		counts <- counts[countsExons$ix,]
	}
	## filter 0 read bins ##
	#counts <- counts[counts$reads > 0, ]
	return(counts)
}

loadHaplotypeAlleleCounts <- function(inCounts, chrs = c(1:22, "X"), fun = "mean",
      minDepth = 10, minNormQual = 200, haplotypeBinSize = 1e5, minSNPsInBin = 3, 
      genomeStyle = "NCBI", sep = "\t", header = TRUE, seqinfo = NULL) {
	if (is.character(inCounts)){
    #### LOAD INPUT READ COUNT DATA ####
    	message("titan: Loading data and phasing information ", inCounts)
    	data <- read.delim(inCounts, header = header, stringsAsFactors = FALSE, 
        		sep = sep)
      colnames(data) <- c("chr", "posn", "refBase", "ref", "nonRefBase", "nonRef", "normQual", "genotype", "phaseSet")
      if (typeof(data[,"posn"])!="integer" || typeof(data[,"ref"])!="integer" || 
          typeof(data[,"nonRef"])!="integer" || is.null(data$genotype) || is.null(data$phaseSet)){
        stop("loadHaplotypeAlleleCounts: Input counts file format does not 
          match required specifications.")		
      }
  }else if (is.data.frame(inCounts)){  #inCounts is a data.frame
    data <- inCounts
  }else{
    stop("loadHaplotypeAlleleCounts: Must provide a filename or data.frame 
      to inCounts")
  }

  # convert to desired genomeStyle and only include autosomes, sex chromosomes
  data[, 1] <- setGenomeStyle(data[, 1], genomeStyle)
   
  ## sort chromosomes
	indChr <- orderSeqlevels(as.character(data[, "chr"]), X.is.sexchrom = TRUE)
	data <- data[indChr, ]
	## sort positions within each chr
	for (x in unique(data[, "chr"])){
		ind <- which(data[, "chr"] == x)
		data[ind, ] <- data[ind[order(data[ind, "posn"], decreasing = FALSE)], ]
	}
  
  ## filter data ##
  data <- cbind(data, start = data$posn, end = data$posn, depth = data$ref + data$nonRef)
  data <- data[data$depth >= minDepth & data$normQual >= minNormQual, ]
  # get max of allele counts
  #data$refOriginal <- data$ref
  #data$ref <- pmax(data$refOriginal, data$nonRef)
  ## use GRanges to find overlap of tiling haplotypeBinSize ##
  phasedAlleles <- getPhasedAlleleFraction(data)
  data$phasedAlleleFraction <- phasedAlleles$allele.fraction
  data$phasedCount <- phasedAlleles$phasedCount
  data$ref.symmetric <- pmax(data$ref, data$nonRef)
  data.gr <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE, seqinfo = seqinfo, ignore.strand = TRUE)
  tile.gr <- unlist(tileGenome(seqinfo, tilewidth = haplotypeBinSize))
  data.gr <- keepSeqlevels(data.gr, chrs)
  tile.gr <- keepSeqlevels(tile.gr, chrs)
  hits <- findOverlaps(query = data.gr, subject = tile.gr)
  data.gr$haplotypeBin <- subjectHits(hits)
  ## use data.table to process haplotype counts by blocks
  data <- as(data.gr, "data.frame")
  data.dt <- data.table(as(data.gr, "data.frame"))
  haploBinSummary <- data.dt[, list(HaplotypeFraction = mean(phasedAlleleFraction), 
      HaplotypeDepth.sum = sum(phasedCount), HaplotypeBinDepth.sum = sum(depth),
      HaplotypeDepth.mean = round(mean(phasedCount)), #round?
      #HaplotypeDepth.mean = round(mean(ref.symmetric)), 
      HaplotypeBinDepth.mean = round(mean(depth)), #round?
      SNPs = length(phasedAlleleFraction)), by = c("phaseSet", "haplotypeBin")]
  # filter bins by number of SNPs #
  haploBinSummary <- haploBinSummary[SNPs >= minSNPsInBin]
  # summary bins with multiple phaseset ID such that haplotypeBin is unique
  haploBinSummary.unique <- haploBinSummary[, list(SNPs = sum(SNPs),
    HaplotypeFraction = sum(HaplotypeFraction * SNPs) / sum(SNPs), 
    HaplotypeDepth.sum = sum(HaplotypeDepth.sum),
    #HaplotypeDepth.sum.symmetric = sum(HaplotypeDepth.sum.symmetric), 
    HaplotypeBinDepth.sum = sum(HaplotypeBinDepth.sum),
    HaplotypeDepth.mean = round(sum(HaplotypeDepth.mean * SNPs) / sum(SNPs)), #round?
    #HaplotypeDepth.mean.symmetric = round(sum(HaplotypeDepth.mean.symmetric * SNPs) / sum(SNPs)), 
    HaplotypeBinDepth.mean = round(sum(HaplotypeBinDepth.mean * SNPs) / sum(SNPs)), #round?
    phaseSet = phaseSet[which.max(SNPs)]), by = haplotypeBin]
    # get symmetric haplotype fraction
  haploBinSummary.unique[, HaplotypeFraction.symmetric := pmax(HaplotypeFraction, 1 - HaplotypeFraction)]
  haploBinSummary.unique[, HaplotypeDepth.sum.symmetric := pmax(HaplotypeDepth.sum, HaplotypeBinDepth.sum - HaplotypeDepth.sum)]
  haploBinSummary.unique[, HaplotypeDepth.mean.symmetric := pmax(HaplotypeDepth.mean, HaplotypeBinDepth.mean - HaplotypeDepth.mean)]
   # set bin column as key so that we can map back to original data
  setkey(haploBinSummary.unique, haplotypeBin) 
  # add the bin summarized values back to data.dt
  data.dt <- cbind(data.dt, select(haploBinSummary.unique[.(subjectHits(hits))], 
      haplotypeBin.aggr = haplotypeBin, HaplotypeFraction.symmetric, HaplotypeDepth.sum,
      HaplotypeDepth.sum.symmetric, HaplotypeBinDepth.sum, HaplotypeDepth.mean, 
      HaplotypeDepth.mean.symmetric, HaplotypeBinDepth.mean, 
      SNPs, phaseSet.aggr = phaseSet))
  data.dt <- na.omit(data.dt)
  data.dt[, phasedCount.haploSymmetric := {
    if (HaplotypeDepth.sum != HaplotypeDepth.sum.symmetric){
      depth - phasedCount 
    }else{
      phasedCount
    }
  }, by=1:nrow(data.dt)]
  alleleData <- select(data.dt, chr=seqnames, posn=start, 
      refOriginal=ref, nonRef=nonRef, tumDepth=depth)
  alleleData$chr <- as.character(alleleData$chr)
  alleleData$ref = pmax(alleleData$refOriginal, alleleData$nonRef)
  haplotypeData <- select(data.dt, chr=seqnames, posn=start, 
      phaseSet=phaseSet.aggr, refOriginal=ref, tumDepthOriginal = depth)
  haplotypeData$chr <- as.character(haplotypeData$chr)
  if (fun == "sum"){
    haplotypeData$ref <- data.dt[, HaplotypeDepth.sum.symmetric]
    haplotypeData$tumDepth <- data.dt[, HaplotypeBinDepth.sum]
    haplotypeData[, HaplotypeRatio := ref / tumDepth]
    #haplotypeData$haplotypeCount <- data.dt[, phasedCount]#data.dt[, HaplotypeDepth.sum]
    haplotypeData$haplotypeCount <- data.dt[, phasedCount.haploSymmetric]
  }else if (fun == "mean"){
    haplotypeData$ref <- data.dt[, HaplotypeDepth.mean.symmetric]
    haplotypeData$tumDepth <- data.dt[, HaplotypeBinDepth.mean]
    haplotypeData[, HaplotypeRatio := ref / tumDepth]
    #haplotypeData$haplotypeCount <- data.dt[, phasedCount]#data.dt[, HaplotypeDepth.mean]
    haplotypeData$haplotypeCount <- data.dt[, phasedCount.haploSymmetric]
  }else if (fun == "SNP"){
    haplotypeData$ref <- data.dt[, phasedCount.haploSymmetric]
    haplotypeData$tumDepth <- data.dt[, depth]
    haplotypeData$HaplotypeRatio <- data.dt[, HaplotypeDepth.sum.symmetric] / data.dt[, HaplotypeBinDepth.sum]
    haplotypeData$haplotypeCount <- data.dt[, phasedCount.haploSymmetric]
  }
  haplotypeData$nonRef <- haplotypeData$tumDepth - haplotypeData$ref
  return(list(haplotypeData=haplotypeData, alleleData=alleleData, data=data.dt))
}

getPhasedAlleleFraction <- function(x, phasedOnly = FALSE){
  gtSplit <- "\\||\\/"  # both 0|1, 1|0, 0/1 (unphased)
  if (phasedOnly){
    gtSplit <- "\\|" # only 0|1, 1|0
  }
  phasedInd <- which(grepl(gtSplit, x$GT))
  altInd <- as.logical(as.numeric(tstrsplit(x$GT[phasedInd], gtSplit)[[1]]))
  refInd <- !altInd
  allele.fraction <- rep(NA, nrow(x))
  allele.fraction[phasedInd[refInd]] <- x[phasedInd[refInd], REF_count] / x[phasedInd[refInd], DEPTH]
  allele.fraction[phasedInd[altInd]] <- x[phasedInd[altInd], ALT_count] / x[phasedInd[altInd], DEPTH]
  phasedCount <- rep(NA, nrow(x))
  phasedCount[phasedInd[refInd]] <- x[phasedInd[refInd], REF_count]
  phasedCount[phasedInd[altInd]] <- x[phasedInd[altInd], ALT_count]
  return(list(allele.fraction = allele.fraction, phasedCount = phasedCount))
}

getHaplotypesFromVCF <- function(vcfFile, chrs = c(1:22, "X"), build = "hg19",
				filterFlags = c("PASS", "10X_RESCUED_MOLECULE_HIGH_DIVERSITY"), 
				minQUAL = 100, minDepth = 10, minVAF = 0.25, altCountField = "AD",
				keepGenotypes = c("1|0", "0|1", "0/1"), snpDB = NULL){
	require(VariantAnnotation)
	#require(data.table)
	message("Loading ", vcfFile)
	vcf <- readVcf(vcfFile, genome = build)
	chrName <- mapSeqlevels(seqlevels(vcf), style="NCBI")
	rowRanges(vcf) <- renameSeqlevels(rowRanges(vcf), na.omit(chrName))
	#keepGenotypes = c("1|0", "0|1", "0/1")
	
	## filter vcf ##
	message("Filtering VCF ...")
	# keep by filter flags
	indFILTER <- rowRanges(vcf)$FILTER %in% filterFlags
	# keep SNPs - ref and alt have length of 1 and only a single allele for ref/alt
	indSNP <- nchar(unstrsplit(CharacterList(rowRanges(vcf)$ALT), sep=",")) == 1 & 
						nchar(unstrsplit(rowRanges(vcf)$REF, sep=",")) == 1
	message("  by quality >=", minQUAL, ")")
	indQUAL <- rowRanges(vcf)$QUAL >= minQUAL
	# keep genotypes
	message("  by genotypes: ", paste0(keepGenotypes, collapse=";"))
	indGeno <- geno(vcf)$GT %in% keepGenotypes
	ind <- indFILTER & indSNP & indQUAL & indGeno
	vcf <- vcf[which(ind)]
	
	message("  by depth (>=", minDepth, ")")
	depth <- geno(vcf)$DP
	indDP <- depth >= minDepth
	message("  by VAF (>=", minVAF, ")")
	minVAF <- min(minVAF, 1 - minVAF) # symmetric from 0-0.5
	if (!is.null(geno(vcf)[[altCountField]])){
	  altCounts <- unlist(lapply(geno(vcf)[[altCountField]], min))
	}else{
	  stop("Alternate read counts not in AD or AO fields.")
	}
	indVAR <- (altCounts / depth) >= minVAF
	vcf <- vcf[which(indDP & indVAR)]
	rm(depth)
	
  if (!is.null(snpDB)){
    message (" by SNP VCF file ", snpDB)
    snp <- readVcf(snpDB, genome = build)
    snpInd <- nchar(unstrsplit(CharacterList(rowRanges(snp)$ALT), sep=",")) == 1 & 
      nchar(unstrsplit(rowRanges(snp)$REF, sep=",")) == 1
    snp <- snp[which(snpInd)]
    hits <- findOverlaps(query = rowRanges(vcf), subject = rowRanges(snp))
    indSNPDB <- queryHits(hits)
    vcf <- vcf[indSNPDB]
  }

	geno.gr <- rowRanges(vcf)	
	values(geno.gr) <- cbind(values(geno.gr), 
			DataFrame(GT = as.character(geno(vcf)$GT), PS = as.character(geno(vcf)$PS)))
	HT <- getPhasedAllele(geno.gr)
	values(geno.gr) <- cbind(values(geno.gr), DataFrame(HT1 = HT$h1, HT2 = HT$h2))
	geno.gr <- keepSeqlevels(geno.gr, chrs)
	
	return(list(vcf.filtered = vcf, geno = geno.gr))
}

# x GRanges object with GT, REF, ALT metadata columns
# returns reference or alternate allele in haplotype
getPhasedAllele <- function(x){
	h1 <- rep(NA, length(x))
	h2 <- rep(NA, length(x))
	# haplotype 1 (h1) ref allele for 1|0
	h1[x$GT == "1|0"] <- as.character(x$REF)[x$GT == "1|0"]
	h2[x$GT == "1|0"] <- as.character(unlist(x$ALT))[x$GT == "1|0"]
	# haplotype 1 (h1) alt allele for 0|1
	h1[x$GT == "0|1"] <- as.character(unlist(x$ALT))[x$GT == "0|1"]	
	h2[x$GT == "0|1"] <- as.character(x$REF)[x$GT == "0|1"]
	return(list(h1 = h1, h2 = h2))
}

plotHaplotypeFraction <- function(dataIn, chr = NULL, type = "HaplotypeRatio", geneAnnot = NULL, 
    spacing = 4,  xlim = NULL, plotType = "haplotypes", ...) {
    if (!type %in% c("HaplotypeRatio", "AllelicRatio")){
      stop("plotHaplotypeFraction: type must be one of 'HaplotypeRatio' or 'AllelicRatio'.")
    }
    if (!plotType %in% c("haplotypes", "copynumber")){
      stop("plotHaplotypeFraction: plotType must be one of 'haplotypes' or 'copynumber'")
    }
    lohCol.hap <- c(`0`="red", `1`="blue")
    lohCol.titan <- c("#00FF00", "#006400", "#0000FF", "#8B0000", 
        "#006400", "#BEBEBE", "#FF0000", "#BEBEBE", 
        "#FF0000")
    names(lohCol.titan) <- c("HOMD", "DLOH", "NLOH", "GAIN", 
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")
    dataIn <- copy(dataIn)
    colnames(dataIn)[1:2] <- c("Chr", "Position")
    #dataIn$AllelicRatio <- pmax(dataIn$refOriginal, dataIn$tumDepthOriginal - dataIn$refOriginal) / dataIn$tumDepthOriginal
    ## set alternating 0 and 1 for phaseSets
    psRuns <- rle(dataIn$PhaseSet)
    ps.id <- 1:length(psRuns$lengths) %% 2
    dataIn$phaseSet.id <- rep(ps.id, psRuns$lengths)
    #dataIn$HaplotypeCount[as.logical(dataIn$phaseSet.id)] <- dataIn$HaplotypeDepth[as.logical(dataIn$phaseSet.id)] - dataIn$HaplotypeCount[as.logical(dataIn$phaseSet.id)]
    #dataIn$HaplotypeRatio[as.logical(dataIn$phaseSet.id)] <- 1 - dataIn$HaplotypeRatio[as.logical(dataIn$phaseSet.id)]
    #dataIn$AllelicRatio <- dataIn$refOriginal / dataIn$tumDepthOriginal
    dataIn[, HaplotypeRatio.1 := HaplotypeRatio]#dataIn$HaplotypeCount / dataIn$HaplotypeDepth
    dataIn[, HaplotypeRatio.2 := 1 - HaplotypeRatio]#(dataIn$HaplotypeDepth - dataIn$HaplotypeCount) / dataIn$HaplotypeDepth
    
    if (!is.null(chr)) {
        for (i in chr) {
            dataByChr <- dataIn[Chr == i, ]
            #dataByChr <- dataByChr[dataByChr[, "TITANcall"] != "OUT", ]
            # plot the data if (outfile!=''){
            # pdf(outfile,width=10,height=6) }
            if (plotType == "haplotypes" || !"TITANcall" %in% colnames(dataByChr)){
              colors.1 <- lohCol.hap[as.character(dataByChr$phaseSet.id)]
              colors.2 <- lohCol.hap[as.character(as.numeric(!dataByChr$phaseSet.id))]
            }else{
              colors.1 <- lohCol.titan[dataByChr[, "TITANcall"]]
              colors.2 <- colors.1
            }
            
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (missing(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), Position]))
            }
            if (type == "HaplotypeRatio"){
              plot(dataByChr[, Position], dataByChr[, HaplotypeRatio.1], 
                  col = colors.1, 
                  pch = 16, xaxt = "n", las = 1, ylab = "Haplotype Fraction", xlim = xlim, 
                  ...)
              points(dataByChr[, Position], dataByChr[, HaplotypeRatio.2], col = colors.2, pch=16, ...)
            }else if (type == "AllelicRatio"){
               plot(dataByChr[, Position], dataByChr[, AllelicRatio], 
                  col = colors.1, 
                  pch = 16, xaxt = "n", las = 1, ylab = "Allelic Fraction", xlim = xlim, 
                  ...)
            }else{
              stop("Need to specify \"type\": HaplotypeRatio or AllelicRatio")
            }
            lines(as.numeric(c(1, dataByChr[nrow(dataByChr), Position])), rep(0.5, 2), type = "l", 
                  col = "grey", lwd = 3)
            
            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        if (plotType == "haplotypes" || !"TITANcall" %in% colnames(dataIn)){
              colors.1 <- lohCol.hap[as.character(dataIn$phaseSet.id)]
              colors.2 <- lohCol.hap[as.character(as.numeric(!dataIn$phaseSet.id))]
            }else{
              colors.1 <- lohCol.titan[dataIn[, TITANcall]]
              colors.2 <- colors.1
        }
        
        # plot for all chromosomes
        coord <- getGenomeWidePositions(dataIn[, Chr], dataIn[, Position])
        if (type == "HaplotypeRatio"){
          plot(coord$posns, as.numeric(dataIn[, HaplotypeRatio.1]), 
            col = colors.1, pch = 16, 
            xaxt = "n", bty = "n", las = 1, ylab = "Haplotype Fraction", ...)
          points(coord$posns, dataIn[, HaplotypeRatio.2], col = colors.1, pch=16, ...)
        }else if (type == "AllelicRatio"){
          plot(coord$posns, as.numeric(dataIn[, AllelicRatio]), 
            col = colors.2, pch = 16, 
            xaxt = "n", bty = "n", las = 1, ylab = "Allelic Fraction", ...)
        }else{
              stop("Need to specify \"type\": HaplotypeRatio or AllelicRatio")
        }
        lines(as.numeric(c(1, coord$posns[length(coord$posns)])), 
            rep(0.5, 2), type = "l", col = "grey", 
            lwd = 3)
        plotChrLines(unique(dataIn[, Chr]), coord$chrBkpt, 
            c(-0.1, 1.1))
        
    }
}

getBXoverlap <- function(sv, bamfile, windowSize = 5000, minReadOverlapSupport = 1, flags = NULL, fields = NULL, tags = c("BX"), bxColName = "bxtags", return.bxtags = FALSE){

  if (is.null(flags)){
    flags <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
          hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
          isFirstMateRead = NA, isSecondMateRead = NA, 
          isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
          isDuplicate = FALSE)
  }
  if (is.null(fields)){
    fields <- scanBamWhat()
  }

  if (length(sv$mateID) == 0){
    sv[, mateID := 1:nrow(sv)]
  }
  overlap.count <- as.data.frame(matrix(NA, nrow = nrow(sv), ncol = 7, dimnames = list(c(sv$mateID), c("Region1","Region1.BXcount","Region2","Region2.BXcount","Support.BXCount","OverlapCount.bxol.only", "OverlapCount.bxol.plus.support"))))
  counts1 <- list()
  counts2 <- list()
  supportBX <- list()
  bxol.plus.support <- list(); bx1.plus.support <- list(); bx2.plus.support <- list()
  for (i in 1:nrow(sv)){
  
    ## region 1 ##
    if (sv$orient_1[i] == "fwd"){
      start1 <- sv$start_1[i] #max(sv$start_1[i] - windowSize/2, 1)
      end1 <- sv$start_1[i] + windowSize
    }else if (sv$orient_1[i] == "rev"){
      end1 <- sv$start_1[i] #sv$start_1[i] + windowSize/2
      start1 <- max(sv$start_1[i] - windowSize, 1)
    }
    region1 <- paste0(sv$chromosome_1[i], ":", start1 , "-", end1)
  
    ## region 2 ##
    if (sv$orient_2[i] == "fwd"){
      start2 <- sv$start_2[i] #max(sv$start_2[i] - windowSize/2, 1)
      end2 <- sv$start_2[i] + windowSize
    }else if (sv$orient_2[i] == "rev"){
      end2 <- sv$start_2[i] #sv$start_2[i] + windowSize/2
      start2 <- max(sv$start_2[i] - windowSize, 1)
    }
    region2 <- paste0(sv$chromosome_2[i], ":", start2, "-", end2)
  
    ## run Rsamtools scanBam for both regions
    regions <- GRanges(c(sv$chromosome_1[i], sv$chromosome_2[i]), 
                      IRanges(c(start1, start2), c(end1, end2)))
    param <- ScanBamParam(flag = flags, simpleCigar = FALSE,
            reverseComplement = FALSE, tag = tags, tagFilter = list(),
            what=fields, which=regions)
    reads <- scanBam(bamfile, param=param)

    ## 1st region
    if (!is.null(reads[[1]]$tag$BX)){
      regionGR1 <- convertBamListToGRanges(reads[[1]])
      regionGR1.contained.pairs <- excludePartialMateMapping(regionGR1, region = regions[1], minMapQ = minMAPQ)
      counts1[[i]] <- as.data.frame(table(regionGR1.contained.pairs$BX), stringsAsFactors=F)
    }else{
      counts1[[i]] <- data.frame(0)
    }
    ## 2nd region
    if (!is.null(reads[[2]]$tag$BX)){
      regionGR2 <- convertBamListToGRanges(reads[[2]])
      regionGR2.contained.pairs <- excludePartialMateMapping(regionGR2, regions[2])
      counts2[[i]] <- as.data.frame(table(regionGR2.contained.pairs$BX), stringsAsFactors=F)
    }else{
      counts2[[i]] <- data.frame(0)
    }
  
    ## find overlap BX ##
    if (i %% 1000 == 0){
      message(i, " breakpoints analyzed")		
    }
    ## check if "bxtags" column exists, then extract BX tags
    #if (sum(!is.na(sv[i,get(bxColName)])) > 0){
    if (length(sv[[bxColName]]) > 0){
      supportBX[[i]] <- sv[i, tstrsplit(tstrsplit(get(bxColName), ","), "_")]$V1
    }else{
      supportBX[[i]] <- NA
    }
    overlap.count[sv$mateID[i], "Support.BXCount"] <- sum(!is.na(supportBX[[i]]))
    #message("Counting BX in ", region1, " and ", region2)
    #overlap.matrix[i,j] <- length(intersect(counts1[counts1[, 2] >= minReadOverlapSupport, 1], counts2[counts2[, 2] >= minReads, 1]))
  	overlap.count[sv$mateID[i], "Region1"] <- region1
    overlap.count[sv$mateID[i], "Region2"] <- region2
    if (nrow(counts1[[i]]) > 1){
      bx1 <- counts1[[i]][counts1[[i]][, 2] >= minReadOverlapSupport, 1]
      bx1.plus.support[[i]] <- union(na.omit(supportBX[[i]]), bx1)
      overlap.count[sv$mateID[i], "Region1.BXcount"] <- length(bx1.plus.support[[i]])
    }else{
      bx1.plus.support[[i]] <- NA
      overlap.count[sv$mateID[i], "Region1.BXcount"] <- 0
    }
    if (nrow(counts2[[i]]) > 1){
      bx2 <- counts2[[i]][counts2[[i]][, 2] >= minReadOverlapSupport, 1]
      bx2.plus.support[[i]] <- union(na.omit(supportBX[[i]]), bx2)
      overlap.count[sv$mateID[i], "Region2.BXcount"] <- length(bx2.plus.support[[i]])
    }else{
      bx2.plus.support[[i]] <- NA
      overlap.count[sv$mateID[i], "Region2.BXcount"] <- 0
    }
    if (nrow(counts1[[i]]) > 1 & nrow(counts2[[i]]) > 1){
      bxol <- intersect(counts1[[i]][counts1[[i]][, 2] >= minReadOverlapSupport, 1], counts2[[i]][counts2[[i]][, 2] >= minReadOverlapSupport, 1])
      bxol.plus.support[[i]] <- union(bxol, na.omit(supportBX[[i]]))
      overlap.count[sv$mateID[i], "OverlapCount.bxol.only"] <- length(bxol)
      overlap.count[sv$mateID[i], "OverlapCount.bxol.plus.support"] <- length(bxol.plus.support[[i]])
    }else{
      bxol.plus.support[[i]] <- NA
      overlap.count[sv$mateID[i], "OverlapCount.bxol.only"] <-  0
      overlap.count[sv$mateID[i], "OverlapCount.bxol.plus.support"] <-  0
    }
  }
  overlap.count <- cbind(mateID = row.names(overlap.count), overlap.count)
  if (return.bxtags){
    return(list(count=as.data.table(overlap.count), bx1=bx1.plus.support, bx2=bx2.plus.support, support=supportBX, overlap=bxol.plus.support))
  }else{
    return(as.data.table(overlap.count))
  }
}

## function to remove overlapped SV from SVABA, LONGRANGER, and GROCSVS.
## removes redundant SV (from other callers other than "tool"
## returns original sv object minus those redundant events.
keepUniqSVcall <- function(sv, tool){
	overlapIdCol <- paste0("overlap.", tool, ".id")
	dupInd <- which(sv[, !is.na(get(overlapIdCol)) & duplicated(get(overlapIdCol), fromLast=T)])
	uniqInd <- setdiff(which(sv[, Tool == tool]), sv[dupInd, get(overlapIdCol)])
	
	rmInd <- c()
	for (i in 1:length(dupInd)){
		toolID <- sv[dupInd[i], get(overlapIdCol)]
		rmInd <- c(rmInd, which(sv[, get(overlapIdCol) == toolID & Tool != tool]))
	}
	keepInd <- setdiff(1:nrow(sv), rmInd)
	return(sv[keepInd])
}


##############################################
############# HELPER FUNCTION ################
##############################################
plotChrLines <- function(chrs, chrBkpt, yrange, cex = 1) {
    # plot vertical chromosome lines
    for (j in 1:length(chrBkpt)) {
        lines(rep(chrBkpt[j], 2), yrange, type = "l", 
            lty = 2, col = "black", lwd = 0.75)
    }
    numLines <- length(chrBkpt)
    mid <- (chrBkpt[1:(numLines - 1)] + chrBkpt[2:numLines])/2
    #chrs[chrs == "X"] <- 23
    #chrs[chrs == "Y"] <- 24
    #chrsToShow <- sort(unique(as.numeric(chrs)))
    #chrsToShow[chrsToShow == 23] <- "X"
    #chrsToShow[chrsToShow == 24] <- "Y"
    axis(side = 1, at = mid, labels = chrs, cex.axis = cex, tick = FALSE)
}

##############################################
############# HELPER FUNCTION ################
##############################################
getGenomeWidePositions <- function(chrs, posns, seqinfo) {
    # create genome coordinate scaffold

    positions <- as.numeric(posns)
    chrsNum <- unique(chrs)
    chrBkpt <- rep(0, length(chrsNum) + 1)
    prevChrPos <- 0
    for (i in 2:length(chrsNum)) {
        chrInd <- which(chrs == chrsNum[i])
        prevChrPos <- seqinfo[i-1, "seqlengths"] + prevChrPos
        #prevChrPos <- positions[chrInd[1] - 1]
        chrBkpt[i] = prevChrPos
        positions[chrInd] = positions[chrInd] + prevChrPos
    }
    chrBkpt[i + 1] <- seqinfo[i, "seqlengths"] + prevChrPos #positions[length(positions)]
    return(list(posns = positions, chrBkpt = chrBkpt))
} 

##############################################
##### FUNCTION TO PLOT GENE FREQUENCY ########
##############################################
## annot must have at least 3 columns: chr, start, stop
## annot must correspond to same rows as lossFreq and gainFreq
plotFrequency <- function(freq, annot, chrs = c(1:22,"X"), seqinfo, ...){
	if (nrow(annot) != length(freq)){
		stop("annot does not have the same number of rows as freq.")
	}
	dots <- list(...)
	
	indChrs <- annot[,1] %in% chrs
	
	## get plotting coords (X-axis) ##
	coord <- getGenomeWidePositions(annot[indChrs,1], annot[indChrs,3], seqinfo)
	
	# setup plot #
	if (is.null(dots$xlim)){
		dots$xlim <- as.numeric(c(1, tail(coord$posns, 1)))
	}
	if (is.null(dots$ylim)){
		dots$ylim <- c(0,max(freq,na.rm=T) + 1)
	}
	if (is.null(dots$col)){
		dots$col <- "red"
	}
	if (is.null(dots$lwd)){
		dots$lwd <- 1
	}
	
	plot(0, type="n", bty="n", xaxt="n", las=2, xlim=dots$xlim, ylim=dots$ylim, ...)

	# plot  frequency #
	lines(coord$posn, freq[indChrs], type = "h", col=dots$col, lwd=dots$lwd)
	indMax <- which(freq==max(freq[indChrs],na.rm=T))
	lines(coord$posn[indMax], freq[indChrs][indMax], type="h", col="red", lwd=dots$lwd*1.5)
	# plot chr lines 
	plotChrLines(unique(annot[indChrs,1]), coord$chrBkpt, dots$ylim)
 	
 	invisible()     	
}


