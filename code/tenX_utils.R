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

#######################################################
### LONG RANGER FUNCTIONS  #####
#######################################################

#######################################################
### Extract SV to data.frame from LONG RANGER VCF #####
#######################################################
getSVfromCollapsedVCF.LR <- function(vcf, chrs = c(1:22, "X"), genomeStyle = "NCBI"){
	if (class(vcf) == "CollapsedVCF"){
			ind.bnd <- which(info(vcf)$SVTYPE == "BND")
			ind.sv <- which(info(vcf)$SVTYPE != "BND")
			message("Processing ", length(ind.bnd), " BND and ", length(ind.sv), " phased records.")
			svGR <- as.data.frame(rowRanges(vcf))
			svGR$ALT <- unlist(svGR$ALT)
			bkpt2 <- str_match(svGR$ALT, "([0-9XY]+):(\\d+)")[, 2:3]
			mateID <- str_match(unlist(info(vcf)$MATEID), "[a-zA-Z0-9]+_[0-9]+")[, 1]
			ID <- names(rowRanges(vcf))
			sv <- data.frame(chromosome_1=as.vector(svGR$seqnames), start_1=as.numeric(svGR$start), 
																chromosome_2=bkpt2[, 1], start_2=as.numeric(bkpt2[, 2]), REF=svGR$REF,
																mateID=NA, alt_1=svGR$ALT, alt_2=NA, FILTER=rowRanges(vcf)$FILTER,
																QUAL=svGR$QUAL, 
																orient_1=NA, orient_2=NA, stringsAsFactors=FALSE)
			## BND ##
			sv[ind.bnd, "mateID"] <- mateID
			## phased ##
			sv[ind.sv, "mateID"] <- rownames(svGR)[ind.sv]
			sv[ind.sv, "chromosome_2"] <- sv[ind.sv, "chromosome_1"]
			sv[ind.sv, "start_2"] <- sv[ind.sv, "start_1"] + info(vcf)$SVLEN[ind.sv]
					
			rownames(sv) <- ID
			# SPAN - need to comput this (-1 for interchr)
			sv <- cbind(sv, SPAN = abs(sv$start_2 - sv$start_1 + 1))
			sv$SPAN[sv$chromosome_1 != sv$chromosome_2] <- -1
			## get breakpoint 2 alt ##
			sv$alt_2 <- sv[ind.bnd, ][sv$mateID, "alt_1"]
      sv <- as.data.table(sv)
      # EVENT - event id
			#if (sum(rownames(info(header(vcf))) %in% c("SVTYPE")) > 0){
			#	sv <- cbind(sv, SVTYPE = info(vcf)$SVTYPE)
			#}			
			if (sum(rownames(info(header(vcf))) %in% c("HAP_ALLELIC_FRAC")) > 0){
				sv <- cbind(sv, HAP_ALLELIC_FRAC = info(vcf)$HAP_ALLELIC_FRAC)
			}
			if (sum(rownames(info(header(vcf))) %in% c("ALLELIC_FRAC")) > 0){
				sv <- cbind(sv, ALLELIC_FRAC = info(vcf)$ALLELIC_FRAC)
			}
			if (sum(rownames(info(header(vcf))) %in% c("PAIRS")) > 0){
				sv <- cbind(sv, DR = info(vcf)$PAIRS)
			}
			if (sum(rownames(info(header(vcf))) %in% c("SPLIT")) > 0){
				sv <- cbind(sv, SR = info(vcf)$SPLIT)
			}
			if (sum(rownames(info(header(vcf))) %in% c("PVAL")) > 0){
				sv <- cbind(sv, PVAL = info(vcf)$PVAL)
			}			
			if (sum(rownames(info(header(vcf))) %in% c("SOURCE")) > 0){
				sv <- cbind(sv, SOURCE = info(vcf)$SOURCE)
			}
			if (sum(rownames(info(header(vcf))) %in% c("PS")) > 0){
				sv <- cbind(sv, PS = info(vcf)$PS)
			}			
			
			# GT - 0 if absent, 1 if present
			if (sum(rownames(geno(header(vcf))) %in% c("GT")) > 0){
			  sv <- cbind(sv, GT = geno(vcf)$GT)
			}		
			
			## filter by chromosome ##
			seqlevelsStyle(chrs) <- genomeStyle
      		sv$chromosome_1 <- mapSeqlevels(sv$chromosome_1, style = genomeStyle)
      		sv$chromosome_2 <- mapSeqlevels(sv$chromosome_2, style = genomeStyle)
      		ind <- sv$chromosome_1 %in% chrs & sv$chromosome_2 %in% chrs
			sv <- sv[ind,]
			
			## get orientation ##
			## check the other alt to see how the current breakpoint is oriented
			orient1 <- str_match(sv$alt_2, "\\]|\\[") 
			orient2 <- str_match(sv$alt_1, "\\]|\\[")
			sv$orient_1[orient1 == "]"] <- "rev"
			sv$orient_1[orient1 == "["] <- "fwd"
			sv$orient_2[orient2 == "]"] <- "rev"
			sv$orient_2[orient2 == "["] <- "fwd"
			
			## for phased SVs: DEL, DUP:TANDEM, INV ##
			sv[alt_1 == "<DEL>", orient_1 := "rev"]; sv[alt_1 == "<DEL>", orient_2 := "fwd"]
			sv[alt_1 == "<DUP:TANDEM>", orient_1 := "fwd"]; sv[alt_1 == "<DUP:TANDEM>", orient_2 := "rev"]
			sv[alt_1 == "<INV>", orient_1 := "rev"]; sv[alt_1 == "<INV>", orient_2 := "rev"]
			
	}else{
		message("vcf object is not a CollapsedVCF")
		sv <- NULL
	}
	
	#sv <- data.table(sv)
	sv <- sv[!is.na(sv$chromosome_1) & !is.na(sv$chromosome_2), ]
	sv <- sv[sv$chromosome_1 %in% chrs & sv$chromosome_2 %in% chrs, ]
	sv <- sv[order(chromosome_1, start_1)]
	sv <- removeDupSV.LR(sv)
	return(sv)
}

getSVfromBEDPE <- function(bedFile, chrs = c(1:22, "X"), skip=0, genomeStyle = "NCBI"){
	bed <- fread(bedFile, skip=skip)	
	orient <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^ORIENT=(.+)"))[2] }))
	type <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^TYPE=(.+)"))[2] }))
	source <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^SOURCE=(.+)"))[2] }))
	allelefrac <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^ALLELIC_FRAC=(.+)"))[2] }))
	hapfrac <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^HAP_ALLELIC_FRAC=(.+)"))[2] }))
	nsplit <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^NSPLIT=(.+)"))[2] }))
	npair <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^NPAIRS=(.+)"))[2] }))
	ps1 <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^PS1=(.+)"))[2] }))
	ps2 <- unlist(lapply(strsplit(bed$info, ";"), function(x){ na.omit(str_match(x, "^PS2=(.+)"))[2] }))
	
	## assign new start1 and start2
	setnames(bed, c("#chrom1", "start1", "chrom2", "start2", "name", "filters"), 
			c("chromosome_1", "start_1", "chromosome_2", "start_2", "mateID", "FILTER"))
	bed[, start_1 := floor(start_1 + (stop1 - start_1)/2)]
	bed[, start_2 := floor(start_2 + (stop2 - start_2)/2)]
	bed[, c("stop1", "stop2") := NULL]

	## assign allele, hap, type to bed ##
	bed[, HAP_ALLELIC_FRAC := hapfrac]
	bed[, ALLELIC_FRAC := allelefrac]
	bed[, DR := npair]
	bed[, SR := nsplit]
	bed[, PS := ps1]; bed[, PS2 := ps2]
	bed[, SOURCE := source]
	bed[, SVTYPE := type]
	bed[, orient := orient]
	bed[chromosome_1 == chromosome_2, SPAN := start_2 - start_1 + 1]
	bed[chromosome_1 != chromosome_2, SPAN := -1]
	
	## assign orientation to bed ##
	bed[SVTYPE=="DEL", orient_1 := "rev"]; bed[SVTYPE=="DEL", orient_2 := "fwd"]
	bed[SVTYPE=="DUP", orient_1 := "fwd"]; bed[SVTYPE=="DUP", orient_2 := "rev"]
	bed[SVTYPE=="INV", orient_1 := "rev"]; bed[SVTYPE=="INV", orient_2 := "rev"]
	bed[orient != "..", orient_1 := sapply(substr(orient, 1, 1), switch, "+"="fwd", "-"="rev", "."=NA)]
	bed[orient != "..", orient_2 := sapply(substr(orient, 2, 2), switch, "+"="fwd", "-"="rev", "."=NA)]
		
	chrs <- as.character(chrs)
	seqlevelsStyle(chrs) <- genomeStyle
	bed$chromosome_1 <- mapSeqlevels(bed$chromosome_1, style = genomeStyle)
	bed$chromosome_2 <- mapSeqlevels(bed$chromosome_2, style = genomeStyle)
    ind <- bed[chromosome_1 %in% chrs & chromosome_2 %in% chrs & !is.na(chromosome_1) & !is.na(chromosome_2), which=TRUE]
	bed <- bed[ind]
	bed <- bed[order(chromosome_1, start_1)]

	#bed[, info := NULL]
	return(copy(bed))
}

#######################################################
### GROCSVS FUNCTIONS  #####
#######################################################
###########################################
### Extract SV to data.frame from VCF #####
###########################################
getSVfromCollapsedVCF.GROC <- function(vcf, chrs = c(1:22, "X")){
	if (class(vcf) == "CollapsedVCF"){
		svGR <- as.data.frame(rowRanges(vcf))
		svGR$ALT <- unlist(svGR$ALT)
		bkpt2 <- str_match(svGR$ALT, "([chr0-9XY]+):(\\d+)")[, 2:3]
		mateID <- info(vcf)$MATEID
		ID <- names(rowRanges(vcf))
		sv <- data.frame(chromosome_1=as.vector(svGR$seqnames), start_1=as.numeric(svGR$start), 
									chromosome_2=bkpt2[, 1], start_2=as.numeric(bkpt2[, 2]), REF=svGR$REF,
									mateID=mateID, alt_1=svGR$ALT, alt_2=NA, FILTER=rowRanges(vcf)$FILTER,
									orient_1=NA, orient_2=NA, stringsAsFactors=FALSE)
		rownames(sv) <- ID
		# SPAN - need to comput this (-1 for interchr)
		sv <- cbind(sv, SPAN = abs(sv$start_2 - sv$start_1 + 1))
		sv$SPAN[sv$chromosome_1 != sv$chromosome_2] <- -1
		## get breakpoint 2 alt ##
		sv$alt_2 <- sv[sv$mateID, "alt_1"]
		sv <- as.data.table(sv)
 		 # EVENT - event id
		if (sum(rownames(info(header(vcf))) %in% c("EVENT")) > 0){
			sv <- cbind(sv, EVENT = as.numeric(info(vcf)$EVENT))
		}
		if (sum(rownames(info(header(vcf))) %in% c("ASSEMBLED")) > 0){
			sv <- cbind(sv, ASSEMBLED = info(vcf)$ASSEMBLED)
		}
		# GT - 0 if absent, 1 if present
		if (sum(rownames(geno(header(vcf))) %in% c("GT")) > 0){
		  sv <- cbind(sv, GT = geno(vcf)$GT)
		}		
		# BS - Number of barcodes supporting the event
		if (sum(rownames(geno(header(vcf))) %in% c("BS")) > 0){
		  sv <- cbind(sv, BS = geno(vcf)$BS)
		}		
		# BT - Total number of barcodes, calculated as the union of barcodes at each site
		if (sum(rownames(geno(header(vcf))) %in% c("BT")) > 0){
		  sv <- cbind(sv, BT = geno(vcf)$BT)
		}			
		# PR - p-value for the event, calculated by resampling from the number of supporting and total barcodes
		if (sum(rownames(geno(header(vcf))) %in% c("PR")) > 0){
			sv <- cbind(sv, PR = geno(vcf)$PR)
		}
		# H1x - String Number of supporting barcodes assigned to haplotype 1 (left side of breakpoint)
		if (sum(rownames(geno(header(vcf))) %in% c("H1x")) > 0){
			sv <- cbind(sv, H1x = geno(vcf)$H1x)
		}
		# H2x - String Number of supporting barcodes assigned to haplotype 2 (left side of breakpoint)
		if (sum(rownames(geno(header(vcf))) %in% c("H2x")) > 0){
			sv <- cbind(sv, H2x = geno(vcf)$H2x)
		}
		# H1y - String Number of supporting barcodes assigned to haplotype 1 (right side of breakpoint)
		if (sum(rownames(geno(header(vcf))) %in% c("H1y")) > 0){
			sv <- cbind(sv, H1y = geno(vcf)$H1y)
		}
		# H2y - String Number of supporting barcodes assigned to haplotype 2 (right side of breakpoint)
		if (sum(rownames(geno(header(vcf))) %in% c("H2y")) > 0){
			sv <- cbind(sv, H2y = geno(vcf)$H2y)
		}
		
		## get orientation ##
		## check the other alt to see how the current breakpoint is oriented
		orient1 <- str_match(sv$alt_2, "\\]|\\[") 
		orient2 <- str_match(sv$alt_1, "\\]|\\[")
		sv$orient_1[orient1 == "]"] <- "rev"
		sv$orient_1[orient1 == "["] <- "fwd"
		sv$orient_2[orient2 == "]"] <- "rev"
		sv$orient_2[orient2 == "["] <- "fwd"
	}else{
		message("vcf object is not a CollapsedVCF")
		sv <- NULL
	}
	
	#sv <- data.table(sv)
	sv <- sv[!is.na(sv$chromosome_1) & !is.na(sv$chromosome_2), ]
	sv <- sv[sv$chromosome_1 %in% chrs & sv$chromosome_2 %in% chrs, ]
	sv <- removeDupSV.GROC(sv)
	sv <- sv[order(chromosome_1, start_1)]
	return(sv)
}

##########################################
########## FILTER VCF BY SUPPORT #########
##########################################
### load GROCSVS vcf file ###
loadGROCSVSVCFtoDataTable <- function(svFile, tumorId, normId = NULL, chrs = c(1:22, "X"),
    filterFlags = NULL, minBXOL = 3, absentBXOL = 1, pValThreshold = 0.1){
  vcf <- tryCatch({
    message("Loading grocsvs vcf ", svFile)
    readVcf(svFile, genome = "hg19")
  }, error = function(x){ 
    return(NA) 
  })
  if (is.null(vcf)){
    stop("Can't find GROCSVS vcf ", svFile)
  }
  sv <- getSVfromCollapsedVCF.GROC(vcf, chrs = chrs)
  if (nrow(sv) == 0){
    stop("No usable breakpoints.")
  }
  ids <- as.vector(na.omit(str_match(colnames(sv), "GT.(.+)")[,2]))
  id <- tumorId
  if (is.null(normId)){
    normId <- setdiff(ids, id)[1] # get the only remaining id which should be normal
  }
  somaticColId <- paste0("SOMATIC.", id)
  supportColId <- paste0("SUPPORT.", id)
  bxolColId <- paste0("BS.",id)
## filter by BXOL support ##
  sv[[supportColId]] <- NA
  sv[[somaticColId]] <- NA
  # filter by FILTER flag #
  if (is.null(filterFlags)){
    filterFlags <- sv[, unique(FILTER)]
  }
  indFilter <- sv$FILTER %in% filterFlags
  # filter by GT (presence in tumor sample but not in normal sample)
  indSomatic <- as.numeric(sv[[paste0("GT.", id)]]) == 1 & as.numeric(sv[[paste0("GT.", normId)]]) == 0
  # filter by supporting barcode counts #
  indBXOL <- sv[[bxolColId]] >= minBXOL & sv[[paste0("BS.",normId)]] <= absentBXOL
  # filter by haplotype overlaps #
  indHaplotype <- (as.numeric(sv[[paste0("H1x.", id)]]) >= minBXOL & as.numeric(sv[[paste0("H1y.", id)]]) >= minBXOL & 
      as.numeric(sv[[paste0("H2x.", id)]]) <= absentBXOL & as.numeric(sv[[paste0("H2y.", id)]]) <= absentBXOL) |
      (as.numeric(sv[[paste0("H1x.", id)]]) <= absentBXOL & as.numeric(sv[[paste0("H1y.", id)]]) <= absentBXOL & 
      as.numeric(sv[[paste0("H2x.", id)]]) >= minBXOL & as.numeric(sv[[paste0("H2y.", id)]]) >= minBXOL)  
  # filter by p-value #
  #sv$PR <= pValThreshold
  sv[[supportColId]][indBXOL] <- "BXOL"
  sv[[supportColId]][indHaplotype] <- "UNIQUE_HAPLOTYPE"
  sv[[somaticColId]][indSomatic] <- TRUE
 
  # sort breakpoint pairs: interchr - chromosome order; intrachr - coordinate order
  sv <- sortBkptPairOrder(sv)
  return(sv)
}


##############################################
########## FIND OVERLAPPING SVS ##############
##############################################
# returns index of x that overlaps y
getOverlapSV <- function(x, y, buffer.x = 100, buffer.y = 100){
	x <- cbind(SV.id = 1:nrow(x), copy(x))

	x.1 <- copy(x)
	setnames(x.1, c("chromosome_1"), c("chr"))
	x.1[, start := start_1 - buffer.x]
	x.1[, end := start_1 + buffer.x]
	x.1 <- as(x.1, "GRanges")
	x.2 <- copy(x)
	setnames(x.2, c("chromosome_2"), c("chr"))
	x.2[, start := start_2 - buffer.x]
	x.2[, end := start_2 + buffer.x]
	x.2 <- as(x.2, "GRanges")
	
	y.1 <- copy(y)
	setnames(y.1, c("chromosome_1"), c("chr"))
	y.1[, start := start_1 - buffer.y]
	y.1[, end := start_1 + buffer.y]
	y.1 <- as(y.1, "GRanges")
	y.2 <- copy(y)
	setnames(y.2, c("chromosome_2"), c("chr"))
	y.2[, start := start_2 - buffer.y]
	y.2[, end := start_2 + buffer.y]
	y.2 <- as(y.2, "GRanges")

	suppressWarnings(hits1 <- findOverlaps(query = x.1, subject = y.1))
	suppressWarnings(hits2 <- findOverlaps(query = x.2, subject = y.2))
	
	overlapInd <- unique(cbind(c(queryHits(hits1), queryHits(hits2)), c(subjectHits(hits1), subjectHits(hits2))))
	overlap.sv <- cbind(x = x[overlapInd[, 1]], y = y[overlapInd[, 2]])
	overlap.sv[, overlap.dist.1 := abs(x.start_1 - y.start_1)]
	overlap.sv[, overlap.dist.2 := abs(x.start_2 - y.start_2)]
	overlap.sv <- overlap.sv[overlap.dist.1 <= buffer.x & overlap.dist.2 <= buffer.x]
	
	x.SV.id <- overlap.sv[, x.SV.id]
	return(list(overlap.ind = x.SV.id, overlap.sv=overlap.sv))	
}

##########################################
##### remove duplicate breakpoints #######
##########################################
removeDupSV.LR <- function(sv){
	bkpts <- unique(str_extract(sv$mateID, "(\\d+)"))
	keepInd <- NULL
	for (i in 1:length(bkpts)){
	  ind <- which(grepl(paste0("call_", bkpts[i], "\\_|call_",bkpts[i], "$"), sv$mateID)) # _47_ | _47$
	  keepInd <- c(keepInd, ind[order(sv[ind, start_1]) == 1])
	}
	sv <- sv[keepInd, ] ## keep 1st mate in pair 
#	 }
	return(sv)
}

##########################################
##### remove duplicate breakpoints #######
##########################################
removeDupSV.GROC <- function(sv){
	bkpts <- unique(str_extract(sv$mateID, "(\\d+)"))
	keepInd <- NULL
	for (i in 1:length(bkpts)){
	  ind <- which(grepl(paste0("^", bkpts[i], ":"), sv$mateID))
	  keepInd <- c(keepInd, ind[order(sv[ind, start_1]) == 1])
	}
	sv <- sv[keepInd, ] ## keep 1st mate in pair 
#	 }
	return(sv)
}


