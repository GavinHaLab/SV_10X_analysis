#' create_circosPlot.R
#' author: Minjeong Ko
#' Fred Hutchinson Cancer Research Center
#' contact: <mko@fredhutch.org>
#' date:  March 19, 2020
#' description: Generate circos plots of copy number and SV breakpoint arcs.

library(RCircos)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID"),
  make_option(c("--svFile"), type = "character", help = "Combined Svaba and Titan SV file (*svabaTitan.sv.PONfilter.bedpe)."),
  make_option(c("--cnFile"), type = "character", help = "Combined Svaba and Titan CNA file (*svabaTitan.cn.txt)."),  
  make_option(c("--genomeBuild"), type="character", default="hg38", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "UCSC", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--excludeCNoverlapType"), type = "character", default = "c(\"Unknown-ShortSVwithCN\")", help = "Exclude SV class."),
  make_option(c("--excludeTools"), type = "character", default = "c()", help = "Exclude SVs predicted by these tools"),
  make_option(c("--outPlotFile"), type="character", help="Path to output figure file.")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

id <- opt$id
svFile <- opt$svFile
cnFile <- opt$cnFile
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
excludeCNoverlapType <- eval(parse(text = opt$excludeCNoverlapType))
excludeTools <- eval(parse(text = opt$excludeTools))
outPlotFile <- opt$outPlotFile

#Load chromosome cytoband data
if(genomeBuild =="hg38"){
	data(UCSC.HG38.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
}else{
	data(UCSC.HG19.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
}
cyto.info <- as.data.table(cyto.info)

#Setup RCircos core components
RCircos.Set.Core.Components(cyto.info, tracks.inside=10, tracks.outside=0) 

# load SV data and filter
sv.data <- fread(svFile, head=T)
if (!is.null(excludeCNoverlapType)){
	sv.data <- sv.data[!CN_overlap_type %in% excludeCNoverlapType]
}
if (!is.null(excludeTools)){
	sv.data <- sv.data[!Tool %in% excludeTools]
}
# check LongRanger coordinates within cyto.info
# if coord is larger than cyto.info, then set it to the largest in cyto.info
chrEnds <- cyto.info[, max(chromEnd), by = Chromosome]
for (i in sv.data[, unique(chrom1)]){
	chrEndCoord <- chrEnds[Chromosome == i, V1]
	sv.data[chrom1 == i & start1 > chrEndCoord, start1 := chrEndCoord]
	sv.data[chrom1 == i & end1 > chrEndCoord, end1 := chrEndCoord]
	sv.data[chrom2 == i & start2 > chrEndCoord, start2 := chrEndCoord]
	sv.data[chrom2 == i & end2 > chrEndCoord, end2 := chrEndCoord]
}

# load CNV data
cnv.data <- fread(cnFile, head=T)
cnv.data.trim <- cnv.data[, c("Chromosome","Start","End","Corrected_Copy_Number")]


pdf(outPlotFile,onefile=T,width=19,height=19)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Histogram.Plot(cnv.data.trim, data.col=4, track.num=1, side="in")
RCircos.Link.Plot(sv.data, track.num=2, by.chromosome=FALSE)
dev.off()
