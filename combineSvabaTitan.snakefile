configfile: "config/config.yaml"
configfile: "config/samples.yaml"

import glob
def getLRFullPath(base, filename):
  return glob.glob(''.join([base, "/*/outs/", filename]))
  
def getTITANpath(base, id, ext):
  return glob.glob(''.join([base, "results/titan/optimalClusterSolution/", id, "_cluster*", ext]))

CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']

#TUM, CLUST = glob_wildcards("../../TITAN/snakemake/results/titan/optimalClusterSolution/{tum}_cluster1.titan.ichor.cna.txt")
#SEG,CLUST = glob_wildcards(config["titanPath"], "/results/titan/optimalClusterSolution/{tumor}_cluster{clust}.titan.ichor.cna.txt")

rule all:
  input: 
  	expand("results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.txt", tumor=config["pairings"]),
  	expand("results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.cn.txt", tumor=config["pairings"]),
  	expand("results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.bedpe", tumor=config["pairings"]),
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt", tumor=config["pairings"]),
  	#expand(directory("results/plotSVABAandTITAN/{tumor}"), tumor=config["pairings"])
  	expand("results/plotSVABAandTITAN/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}", tumor=config["pairings"], type=config["plot_type"], chr=CHRS, format=config["plot_format"])
  		
rule getLongRangerSomaticSV:
	input:
		tumSVFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "large_sv_calls.bedpe"),
		normSVFile=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], "large_sv_calls.bedpe"),
		tumDelFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "dels.vcf.gz"),
		normDelFile=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], "dels.vcf.gz")
	output:
		outputSVFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt"
	params:
		getLRscript=config["getLRsomaticSV_script"],
		tenXfuncs=config["tenX_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"]
	log:
		"logs/LongRangerSomaticSV/{tumor}.log"
	shell:
		"Rscript {params.getLRscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --tumLargeSVFile {input.tumSVFile} --normLargeSVFile {input.normSVFile} --tumDeletionFile {input.tumDelFile} --normDeletionFile {input.normDelFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/LongRangerSomaticSV/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} > {log} 2> {log}"

rule combineSVABAandTITAN:
	input:
		LRsummaryFile=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], "summary.csv"),
		svabaVCF="results/barcodeRescue/{tumor}.bxOverlap.vcf",
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		LRsvFile="results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt"
	output:
		outputSVFile="results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.txt",
		outputBedpeFile="results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.bedpe",
		outputCNFile="results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.cn.txt"
	params:
		combineSVCNscript=config["combineSVCN_script"],
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
		"logs/combineSVABAandTITAN/{tumor}.log"
	shell:
		"Rscript {params.combineSVCNscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --svabaVCF {input.svabaVCF} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --LRsummaryFile {input.LRsummaryFile} --LRsvFile {input.LRsvFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/combineSVABAandTITAN/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} --outputCNFile {output.outputCNFile} --outputBedpeFile {output.outputBedpeFile} > {log} 2> {log}"
		
rule plotSVABAandTITAN:
	input:
		svabaVCF="results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.txt",
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		titanParamFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".params.txt")
	output:
		"results/plotSVABAandTITAN/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}"
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
		"logs/plotSVABAandTITAN/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}.log"
	shell:
		"Rscript {params.plotSVCNscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --plot_funcs {params.plotfuncs} --titan_libdir {params.libdir} --svFile {input.svabaVCF} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --titanParamFile {input.titanParamFile} --chrs \"{params.chrs}\" --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --cytobandFile {params.cytobandFile} --start {params.start} --end {params.end} --zoom {params.zoom} --plotYlim \"{params.ylim}\" --geneFile {params.geneFile} --plotCNAtype {params.type} --plotSize \"{params.size}\" --outPlotFile {output} > {log} 2> {log}" 
	

