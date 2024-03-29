#' plotTitanSvaba.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:  August 3, 2018
#' description: Generate plots of copy number and SV breakpoint arcs.


library(optparse)
option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID"),
  make_option(c("--tenX_funcs"), type = "character", help = "Path to file containing 10X R functions to source."),
  make_option(c("--svaba_funcs"), type = "character", help = "Path to file containing SVABA R functions to source."),
  make_option(c("--plot_funcs"), type = "character", help = "Path to file containing plotting R functions to source."),
  make_option(c("--titan_libdir"), type = "character", help = "Directory containing source code. Specify if changes have been made to source code and want to over-ride package code."),
  make_option(c("--svFile"), type="character", help = "Combined Svaba and Titan SV file (*svabaTitan.sv.txt)."),
  make_option(c("--titanBinFile"), type="character", help = "Path to TITAN titan.txt output file."),
  make_option(c("--titanSegFile"), type="character", help = "Path to TITAN segs.txt output file."),
  make_option(c("--titanParamFile"), type="character", help = "Path to TITAN params.txt output file."),
  make_option(c("--geneFile"), type="character", default=NULL, help = "Path to file containing list of genes with chr, start, end coordinates."),
  make_option(c("--customSVFile"), type="character", default=NULL, help = "Path to additional SV events (not in --svFile) to plot."),
  make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--cytobandFile"), type = "character", default = NULL, help = "Cytoband file should be provided only if reference genome is hg38."),
  make_option(c("--zoom"), type = "logical", default = FALSE, help = "Zoom plot; if TRUE, then requires --chrs --start --end to be set. [Default: %default]"),
  make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to plot; string [Default: %default"),
  make_option(c("--start"), type = "integer", default = NULL, help = "Start coordinate for zoom plots"),
  make_option(c("--end"), type = "integer", default = NULL, help = "End coordinate for zoom plots"),
  make_option(c("--yaxis"), type = "character", default = "logratio", help = "The y-axis data type: integer or logratio. [Default: %default]."),
  make_option(c("--plotCNAtype"), type="character", default = "titan", help = "titan or ichor; if titan, then will also plot haplotype fraction. [Default: %default]."),
  make_option(c("--plotYlim"), type = "character", default = "c(-2,6)", help = "Y limits for plotting log ratio. [Default: %default]."),
  make_option(c("--plotSize"), type = "character", default = "c(8,4)", help = "width and height in inches. [Default: %default]."),
  #make_option(c("--plotFormat"), type = "character", default = "png", help = "File format of plot. E.g. pdf or png. [Default: %default]."),
  make_option(c("--outPlotFile"), type="character", help="Path to output figure file.")
  #make_option(c("--outDir"), type="character", help="Path to output directory.")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

options(bitmapType='cairo', scipen=0)
options(stringsAsFactors=F, bitmapType = "cairo", width=175)

library(data.table)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(reshape2)
library(diagram)
library(tools)
library(SNPchip)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
setDTthreads(1)
source(paste0(opt$titan_libdir, "/R/haplotype.R"))
#source(opt$tenX_funcs)
source(opt$svaba_funcs)
source(opt$plot_funcs)


id <- opt$id
svFile <- opt$svFile
cnFile <- opt$titanBinFile
segFile <- opt$titanSegFile
paramFile <- opt$titanParamFile
chrStr <- as.character(eval(parse(text = "c(opt$chrs)")))
startPos <- opt$start
endPos <- opt$end
zoom <- opt$zoom
ylim <- eval(parse(text = opt$plotYlim))
geneList <- opt$geneFile
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle 
cytobandFile <- opt$cytobandFile
altSVFile <- opt$customSVFile
plotType <- opt$plotCNAtype
plotSize <- eval(parse(text=opt$plotSize))
#plotFormat <- opt$plotFormat
#outDir <- opt$outDir
outPlotFile <- opt$outPlotFile
plotFormat <- tools::file_ext(outPlotFile)
width <- plotSize[1]  #6 8 
height <- plotSize[2]  #3 3.5 #4 
spacing <- 3
yaxis <- opt$yaxis #"integer"
minSPAN <- 10000
ylimSV <- ylim
plotArrows <- FALSE
plotSegs <- FALSE
interchr <- TRUE
plotHaplotypeFrac <- FALSE
buffer <- 1e4
offset.factor <- 1.15 # sv drawn outside of plot
lcol <- NULL
arr.col <- NULL
svabaCol <- "black"
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
svCol <- c('TandemDuplication'='#D55E00', 'Deletion'='#009E73', 'Inversion'='#332288',
  'Translocation'='#0072B2','Balanced'='#F0E442','Unbalanced'='#CC79A7') #"black"
rescueCol <- "purple"
lrCol <- "lightgrey"
grocCol <- "blue"
manualCol <- "purple"
#seqinfo <- Seqinfo(genome=genomeBuild)
bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
	seqinfo <- Seqinfo(genome=genomeBuild)
} else {
	seqinfo <- seqinfo(get(bsg))
}


if (zoom){
  xlim <- c(startPos, endPos)
  startTitle <- paste0(format(round(startPos/1e6,2), nsmall=2))
  endTitle <- paste0(format(round(endPos/1e6,2),nsmall=2), "Mb")
  cex <- 0.75
  cytoBand <- F
  xaxt <- "s"
  plotAtCentre <- FALSE
  cnColor <- FALSE
  plotIdio <- FALSE
  exclude.na.snp <- FALSE
}else{
  xlim <- NULL
  cex <- 0.25
  cytoBand <- T
  xaxt <- "n"
  #yaxis <- "logratio"
  plotAtCentre <- FALSE
  cnColor <- FALSE
  plotIdio <- TRUE
  plotSegs <- FALSE
  exclude.na.snp <- TRUE
}
if (chrStr == "0" || is.null(chrStr) || chrStr == "None"){
  chrStr <- as.character(c(1:22, "X"))
}
seqlevelsStyle(chrStr) <- genomeStyle

if (plotType == "titan"){
	cnColor <- TRUE
	plotHaplotypeFrac <- TRUE
	height <- height * 1.5
	cex.hap <- 0.5
}

if (!cnColor){
  cnCol <- rep("#000000", 30)
}else{
  cnCol <- NULL
}

if (!is.null(altSVFile) && altSVFile != "None"){
	altSV <- fread(altSVFile, na.string = ".")
}

if (!is.null(geneList) && geneList != "None"){
  genes <- read.delim(geneList, header=F, as.is=T)
}else{
  genes <- NULL
}
colnames(genes) <- c("Gene", "Chr", "Start", "End")
seqlevelsStyle(genes$Chr) <- genomeStyle

if (genomeBuild == "hg38" && file.exists(cytobandFile)){
  cytoband <- as.data.frame(fread(cytobandFile))
  names(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
}

#outPlotDir <- outDir
#dir.create(outPlotDir)
#outImage <- paste0(outDir, "/", id, ".RData")
outImage <- gsub(plotFormat, "RData", outPlotFile)
#save.image(file=outImage)

message("Analyzing ", id)
ulp <- fread(cnFile)
segs <- fread(segFile)
#segs$median <- segs$seg.mean
#segs$chr <- segs$chrom
#segs$state <- segs$state.num
params <- read.delim(grep(id, paramFile, value=T), header=F, as.is=T, nrow=2)
purity <- 1 - as.numeric(params[1, 2]) # TITAN normal contamination
ploidyT <- as.numeric(params[2, 2])
normCN <- 2
ploidyS <- purity * ploidyT + (1-purity) * normCN
if (yaxis == "integer"){
	ulp[!grepl("X",Chr), LogRatio := log2(logRbasedCN(LogRatio, purity, ploidyT, cn=2))]
	ulp[grepl("X",Chr), LogRatio := log2(logRbasedCN(LogRatio, purity, ploidyT, cn=1))]
	colName <- "logR_Copy_Number"
}else{
	ulp[, LogRatio := LogRatio + log2(ploidyS / 2)]
	segs$LogRatio <- segs$Median_logR
	segs$LogRatio <- segs$LogRatio + log2(ploidyS / 2)
	colName <- "LogRatio"
}
# exclude data points not analyzed by titan
if (plotType == "titan"){
	ulp <- ulp[!is.na(Corrected_Call)]
}else{
  ulp[is.na(Corrected_Call), Corrected_Call := "NEUT"]
}
ulp <- ulp[Chr %in% chrStr]
ulp$Chr <- factor(ulp$Chr, levels = chrStr)
ulp <- ulp[order(Chr, Start)]

## for plotting full chromosomes, exclude rows where it is not a SNP
# usually, these rows are ichorCNA bins and can contain noisy data
# for zoom, we leave some of these in to avoid gaps in the plots
if (exclude.na.snp){
	ulp <- ulp[!is.na(Position)]
}

############# load Combined SV (SVABA, GROC, LongRanger) ##############
sv <- fread(svFile, na.string = ".")
# filter very short SVABA calls with no CN_overlap_type or short SVABA calls with CN overlap
sv <- sv[SV.Filter == "KEEP" & !is.na(CN_overlap_type)]
#sv <- sv[CN_overlap_type != "." & CN_overlap_type != "Unknown-ShortSVwithCN"]
# consolidate SV classes into those defined by `svabaCol`
#Tandem Duplication           Deletion          Inversion              Trans           Balanced         Unbalanced
sv[, SV.class := CN_overlap_type]
sv[grepl("TandemDup", CN_overlap_type), SV.class := "TandemDuplication"]
sv[grepl("Inversion", CN_overlap_type), SV.class := "Inversion"]
sv[grepl("Trans", CN_overlap_type), SV.class := "Translocation"]
sv[is.na(CN_overlap_type), SV.class := "Unbalanced"]
sv[, color := svCol[SV.class]]
#save.image(file=outImage)

#####################################
########## PLOT CHR RESULTS #########
#####################################	 
for (j in 1:length(chrStr)){
  message("Plotting ", chrStr[j])
  ###################################
  ### SNOWMAN + BARCODE RESCUE ######
  if (genomeStyle == "NCBI"){
  	chrTitle <- paste0("chr", chrStr[j])
  }else{
  	chrTitle <- chrStr[j]
  }
  
  #outPlot <- paste0(outPlotDir, "/", id, "_CNA-SV-BX_",plotType,"_",chrTitle,".",plotFormat)
  plotTitle <- paste0(id, " (", chrTitle,")")
  if (zoom){
    ylimMax <- ulp[Chr==chrStr[j] & Start >= xlim[1] & Start <= xlim[2], max(get(colName), na.rm=T)] + 1
    #outPlot <- paste0(outPlotDir, "/", id, "_CNA-SV-BX_",plotType,"_",chrTitle,"-",startTitle,"-",endTitle,".",plotFormat)
    plotTitle <- paste0(id, " (",chrTitle,":",startTitle,"-",endTitle, ")")
  }else{
    xlim <- c(1, seqlengths(seqinfo)[chrStr[j]])
    ylimMax <- ulp[, max(get(colName), na.rm=T)] + 1
  }    
  ylim[2] <- min(max(ylim[2], ceiling(ylimMax)), 10)
  ylimSV <- ylim
  ylimSV[2] <- ylimSV[2] - 0.5
  
  if (plotFormat == "png"){
  	png(outPlotFile, width = width*100, height=height*100)
  }else{
  	pdf(outPlotFile, width = width, height=height)
	}
  if (plotHaplotypeFrac){ 
    par(mfrow=c(2,1)); spacing <- 0  
  }else{
    spacing <- 6
  }	
  if (plotSegs) { segsToPlot <- segs } else { segsToPlot <- NULL}
  
  if (grepl("X", chrStr[j])) { cnCol <- rep("#000000", 30) }
  message("Plotting read depth CN")
  plotTitanIchorCNA(ulp, segs=segsToPlot, chr=chrStr[j], colName=colName, 
      cytoBand=FALSE, geneAnnot=genes, purity = purity, ploidyT = NULL, yaxis=yaxis, cnCol = cnCol,
      yrange=ylim, xlim=xlim, spacing=spacing, xaxt=xaxt, cex = cex, gene.cex = 1,
      plot.title = plotTitle)

  if (nrow(sv) > 0){
      centreLine <- 0
    
    if (yaxis == "integer"){
   		normCN <- ifelse(grepl("X", chrStr[j]), 1, 2) 
    	ylimMax.int <- ulp[Chr == chrStr[j], max(get(colName), na.rm=T)] * 0.75
    	ylimSV[2] <- min(ylimMax.int, 10)
    	ylimSV[1] <- 2
    	centreLine <- ifelse(grepl("X", chrStr[j]), 0, 1)
    	ploidyS <- NULL
    }
    #arcHeight <- ylim[2]# - ylim[1])/4
    
    if (!is.null(altSVFile) && altSVFile != "None"){
			altSV.sample <- altSV[Sample == id]
			message("Plotting custom")
			plotRearrangementArcs(altSV.sample, cn=as.data.frame(ulp), chr=chrStr[j], 
										interchr = interchr, plotAtCentre = plotAtCentre,
										xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
										centreLine=centreLine, buffer=buffer, lcol=manualCol, arr.col=manualCol, 
										endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
  	}	

    ## plot LongRanger arcs ##
    message("Plotting longranger")
    plotRearrangementArcs(sv[Tool=="LONGRANGER"], cn=as.data.frame(ulp), chr=chrStr[j], 
                  interchr = interchr, plotAtCentre = plotAtCentre,
                  xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
                  centreLine=centreLine, buffer=buffer, svTypeCol=!is.null(svCol), 
                  lcol=lrCol, arr.col=lrCol, lwd = 2,
                  endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
    ## plot BX Rescue arcs ##
    message("Plotting CN rescue")
    plotRearrangementArcs(sv[Tool=="SVABA" & support %in% c("RESCUE,CN2","RESCUE,CN2") &
                          CN_overlap_type != "Unknown-ShortSVwithCN"], 
    							cn=as.data.frame(ulp), 
    							chr=chrStr[j], interchr = interchr, plotAtCentre = plotAtCentre,
                  xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
                  centreLine=centreLine, buffer=buffer, svTypeCol=!is.null(svCol),
                  lcol=rescueCol, arr.col=rescueCol, lwd = 2,
                  endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
    # message("Plotting CN rescue - short SVs")
    # plotRearrangementArcs(sv[Tool=="SVABA" & support %in% c("BX,CN1","BX,CN2","SVABA,CN1","SVABA,CN2") &
    #                       CN_overlap_type == "Unknown-ShortSVwithCN"], 
    #               cn=as.data.frame(ulp), 
    #               chr=chrStr[j], interchr = interchr, plotAtCentre = plotAtCentre,
    #               xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
    #               centreLine=centreLine, buffer=buffer, lcol=rescueCol, arr.col=rescueCol, 
    #               endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
    message("Plotting BX rescue")
    plotRearrangementArcs(sv[Tool=="SVABA" & (support %in% c("BX,CN1","BX,CN2","BX"))], 
    							cn=as.data.frame(ulp), 
    							chr=chrStr[j], interchr = interchr, plotAtCentre = plotAtCentre,
                  xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
                  centreLine=centreLine, buffer=buffer, svTypeCol=!is.null(svCol), 
                  lcol=rescueCol, arr.col=rescueCol, lwd = 2,
                  endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
    message("Plotting manual SVs")
    plotRearrangementArcs(sv[Tool=="SVABA" & (support %in% c("MANUAL","MANUAL,CN1","MANUAL,CN2"))], 
                  cn=as.data.frame(ulp), 
                  chr=chrStr[j], interchr = interchr, plotAtCentre = plotAtCentre,
                  xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
                  centreLine=centreLine, buffer=buffer, svTypeCol=!is.null(svCol),
                  lcol=manualCol, arr.col=manualCol, lwd = 2,
                  endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
    ## plot SVABA arcs ##
    message("Plotting svaba")
    plotRearrangementArcs(sv[Tool=="SVABA" & support %in% c("SVABA","SVABA,CN1","SVABA,CN2")], cn=as.data.frame(ulp), 
    							chr=chrStr[j], interchr = interchr, plotAtCentre = plotAtCentre,
                  xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
                  centreLine=centreLine, buffer=buffer, svTypeCol=!is.null(svCol),
                  lcol=svabaCol, arr.col=svabaCol, lwd = 2,
                  endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
    ## plot GROCSVS arcs ##
    message("Plotting grocsvs")
    plotRearrangementArcs(sv[Tool=="GROCSVS"], cn=as.data.frame(ulp), chr=chrStr[j], 
                  interchr = interchr, plotAtCentre = plotAtCentre,
                  xlim=xlim, arcHeight=ylimSV, ploidy = NULL, lty = 1, offset.factor=offset.factor,
                  centreLine=centreLine, buffer=buffer, svTypeCol=!is.null(svCol),
                  lcol=grocCol, arr.col=grocCol, lwd = 2,
                  endhead = plotArrows, arr.pos = 1.0, minSPAN = 0)
  }
  
  if (plotHaplotypeFrac){
    message("Plotting haplotype fraction")
    plotHaplotypeFraction(ulp[,-1], chrStr[j], resultType = "HaplotypeRatio", colType = "Haplotypes", 
	  xlab="", ylim=c(0,1), xlim=xlim, cex=cex.hap, cex.axis=1.5, cex.lab=1.5, spacing = 6)
	
	par(xpd=NA)
    
    if (genomeBuild == "hg38" && file.exists(cytobandFile)){
      sl <- seqlengths(seqinfo[chrStr[j]])
      pI <- plotIdiogram.hg38(chrStr[j], cytoband=cytoband, seqinfo=seqinfo, xlim=c(0, max(sl)), unit="bp", label.y=-0.425, new=FALSE, ylim=c(-0.3,-0.15))	
    }else{
      pI <- plotIdiogram(chrStr[j], build="hg19", unit="bp", label.y=-0.6, new=FALSE, ylim=c(-0.3,-0.15))
    }
  }else{ # not plotting haplotype fraction
    if (!zoom){
      par(xpd=NA)
      if (genomeBuild == "hg38" && file.exists(cytobandFile)){
        sl <- seqlengths(seqinfo[chrStr[j]])
        pI <- plotIdiogram.hg38(chrStr[j], cytoband=cytoband, seqinfo=seqinfo, unit="bp", label.y=ylim[1]-(ylim[2]-ylim[1])*0.275, new=FALSE, ylim=c(ylim[1]-(ylim[2]-ylim[1])*0.15,ylim[1]-(ylim[2]-ylim[1])*0.075))
      }else{
        pI <- plotIdiogram(chrStr[j], build="hg19", unit="bp", label.y=ylim[1]-(ylim[2]-ylim[1])*0.275, new=FALSE, ylim=c(ylim[1]-(ylim[2]-ylim[1])*0.15,ylim[1]-(ylim[2]-ylim[1])*0.075))
      }
    }
  }

  
  dev.off()
		
}
