configfile: "config/config.yaml"
configfile: "config/samples.yaml"

import glob
import re
def getLRFullPath(base, filename):
  return glob.glob(''.join([base, "/*/outs/", filename]))
  
def getTITANpath(base, id, ext):
  return glob.glob(''.join([base, "results/titan/optimalClusterSolution/", id, "_cluster*", ext]))

def getGROCpath(base, id):
  m = re.search('[0-9]+', id)
  return glob.glob(''.join([base, m.group(), "/results/PostprocessingStep/svs.vcf"]))

CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']

#TUM, CLUST = glob_wildcards("../../TITAN/snakemake/results/titan/optimalClusterSolution/{tum}_cluster1.titan.ichor.cna.txt")
#SEG,CLUST = glob_wildcards(config["titanPath"], "/results/titan/optimalClusterSolution/{tumor}_cluster{clust}.titan.ichor.cna.txt")

rule all:
  input: 
  	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt", tumor=config["pairings"]),
  	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.cn.txt", tumor=config["pairings"]),
  	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.bedpe", tumor=config["pairings"]),
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt", tumor=config["pairings"]),
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.germline.sv.txt", tumor=config["pairings"]),
  	"results/panelOfNormalsSV/PanelOfNormalsSV.txt",
	"results/panelOfNormalsSV/PoNBlacklistBins.txt",
	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe", tumor=config["pairings"]),
	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNfilter.bedpe", tumor=config["pairings"]),
  	expand("results/plotSvabaGrocsvsTitan/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}", tumor=config["pairings"], type=config["plot_type"], chr=CHRS, format=config["plot_format"])
  		
rule getLongRangerSomaticSV:
	input:
		tumSVFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "large_sv_calls.bedpe"),
		normSVFile=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], "large_sv_calls.bedpe"),
		tumDelFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "dels.vcf.gz"),
		normDelFile=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], "dels.vcf.gz")
	output:
		outputSVFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt",
		outputNormSVFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.germline.sv.txt",
	params:
		getLRscript=config["getLRsomaticSV_script"],		
		tenXfuncs=config["tenX_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"]
	log:
		"logs/LongRangerSomaticSV/{tumor}.log"
	shell:
		"Rscript {params.getLRscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --tumLargeSVFile {input.tumSVFile} --normLargeSVFile {input.normSVFile} --tumDeletionFile {input.tumDelFile} --normDeletionFile {input.normDelFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/LongRangerSomaticSV/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} --outputNormSVFile {output.outputNormSVFile} > {log} 2> {log}"

rule buildPoN:
	input:
		svabaDir="results/svaba/",
		lrDir="results/LongRangerSomaticSV/"
	output:
		outputPoNFile="results/panelOfNormalsSV/PanelOfNormalsSV.txt",
		outputBlackListFile="results/panelOfNormalsSV/PoNBlacklistBins.txt"
	params:
		buildPoNscript=config["buildPoN_script"],
		blackListBinWidth=config["PoN_blackListBinWidth"],
		svabafuncs=config["svaba_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"]
	log:
		"logs/panelOfNormalsSV/panelOfNormalsSV.log"
	shell:
		"Rscript {params.buildPoNscript} --SVABAdir {input.svabaDir} --LRdir {input.lrDir} --svaba_funcs {params.svabafuncs} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outputPoNFile {output.outputPoNFile} --outputBlackListFile {output.outputBlackListFile} > {log} 2> {log}"

rule combineSvabaGrocsvsTitan:
	input:
		LRsummaryFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "summary.csv"),
		svabaVCF="results/barcodeRescue/{tumor}.bxOverlap.vcf",
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		LRsvFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt",
		grocFile=lambda wildcards: getGROCpath(config["grocsvs_results"], wildcards.tumor)
	output:
		outputSVFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt",
		outputBedpeFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.bedpe",
		outputCNFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.cn.txt"
	params:
		combineSVCNscript=config["combineSVCN_script"],
		normID=lambda wildcards: config["pairings"][wildcards.tumor],
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"],
		minMapQ=config["bxRescue_minMapQ"],
		minLength=config["bxRescue_minLength"],
		windowSize=config["bxRescue_windowSize"],
		minRead=config["bxRescue_minReadOverlapSupport"]	
	log:
		"logs/combineSvabaGrocsvsTitan/{tumor}.log"
	shell:
		"Rscript {params.combineSVCNscript} --tumID {wildcards.tumor} --normID {params.normID} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --svabaVCF {input.svabaVCF} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --LRsummaryFile {input.LRsummaryFile} --LRsvFile {input.LRsvFile} --grocsvsFile {input.grocFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/combineSvabaGrocsvsTitan/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} --outputCNFile {output.outputCNFile} --outputBedpeFile {output.outputBedpeFile} > {log} 2> {log}"

rule filterPoNSV:
	input:
		svFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt",
		PoNFile="results/panelOfNormalsSV/PanelOfNormalsSV.txt",
		blackListFile="results/panelOfNormalsSV/PoNBlacklistBins.txt"
	output:
		outputSVAnnotFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe",
		outputSVFiltFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNfilter.bedpe",
		outputSummaryFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNfilter.rmCounts.txt"
	params:
		filterScript=config["filterPoNSV_script"],
		svabafuncs=config["svaba_funcs"],
		minFreqPoNSVBkptOverlap=config["PoN_minFreqSVbkpts"],
		minFreqPoNBlackList=config["PoN_minFreqBlackList"]
	log:
		"logs/combineSvabaGrocsvsTitan/{tumor}.filterPoNSV.log"
	shell:
		"Rscript {params.filterScript} --id {wildcards.tumor} --svaba_funcs {params.svabafuncs} --svFile {input.svFile} --PoNFile {input.PoNFile} --blackListFile {input.blackListFile} --minFreqPoNSVBkptOverlap {params.minFreqPoNSVBkptOverlap} --minFreqPoNBlackList {params.minFreqPoNBlackList} --outputSVAnnotFile {output.outputSVAnnotFile} --outputSVFiltFile {output.outputSVFiltFile} --outputSummary {output.outputSummaryFile} 2> {log} > {log}"

rule plotSvabaGrocsvsTitan:
	input:
		svabaVCF="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNfilter.bedpe",
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		titanParamFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".params.txt")
	output:
		"results/plotSvabaGrocsvsTitan/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}"
	params:
		plotSVCNscript=config["plotSVCN_script"],
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		plotfuncs=config["plot_funcs"],
		libdir=config["titan_libdir"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		cytobandFile=config["cytobandFile"],
		zoom=config["plot_zoom"],
		chrs=config["plot_chrs"],
		start=config["plot_startPos"],
		end=config["plot_endPos"],
		ylim=config["plot_ylim"],
		type=config["plot_type"],
		geneFile=config["plot_geneFile"],
		size=config["plot_size"],
		format=config["plot_format"]	
	log:
		"logs/plotSvabaGrocsvsTitan/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}.log"
	shell:
		"Rscript {params.plotSVCNscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --plot_funcs {params.plotfuncs} --titan_libdir {params.libdir} --svFile {input.svabaVCF} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --titanParamFile {input.titanParamFile} --chrs {wildcards.chr} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --cytobandFile {params.cytobandFile} --start {params.start} --end {params.end} --zoom {params.zoom} --plotYlim \"{params.ylim}\" --geneFile {params.geneFile} --plotCNAtype {params.type} --plotSize \"{params.size}\" --outPlotFile {output} > {log} 2> {log}" 
	

