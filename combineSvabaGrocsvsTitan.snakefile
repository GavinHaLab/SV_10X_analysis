"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Python/3.7.4-foss-2019b-fh1
ml SAMtools/1.10-GCCcore-8.3.0

snakemake -s combineSvabaGrocsvsTitan.snakefile --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50 -np
"""

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
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.somatic.sv.txt", tumor=config["pairings"]),
  	expand("results/LongRangerSomaticSV/{tumor}/{tumor}.LR.germline.sv.txt", tumor=config["pairings"]),
  	"results/panelOfNormalsSV/PanelOfNormalsSV.txt",
	"results/panelOfNormalsSV/PoNBlacklistBins.txt",
  	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt", tumor=config["pairings"]),
  	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.cn.txt", tumor=config["pairings"]),
  	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.bedpe", tumor=config["pairings"]),
 	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe", tumor=config["pairings"]),
	expand("results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.bedpe", tumor=config["pairings"]),
  	expand("results/plotSvabaGrocsvsTitan/{tumor}/{tumor}_CNA-SV-BX_{type}_chr{chr}.{format}", tumor=config["pairings"], type=config["plot_type"], chr=CHRS, format=config["plot_format"]),
   	expand("results/plotCircos/{tumor}/{tumor}_Circos.pdf", tumor=config["pairings"])
 		
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
		grocFile=lambda wildcards: getGROCpath(config["grocsvs_results"], wildcards.tumor),
	output:
		outputSVFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt",
		outputBedpeFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.bedpe",
		outputCNFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.cn.txt"
	params:
		combineSVCNscript=config["combineSVCN_script"],
		normID=lambda wildcards: config["pairings"][wildcards.tumor],
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		manualSVfile=config["manualSVFile"],
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
		"Rscript {params.combineSVCNscript} --tumID {wildcards.tumor} --normID {params.normID} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --svabaVCF {input.svabaVCF} --manualSVFile {params.manualSVfile} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --LRsummaryFile {input.LRsummaryFile} --LRsvFile {input.LRsvFile} --grocsvsFile {input.grocFile} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --outDir results/combineSvabaGrocsvsTitan/{wildcards.tumor}/ --outputSVFile {output.outputSVFile} --outputCNFile {output.outputCNFile} --outputBedpeFile {output.outputBedpeFile} > {log} 2> {log}"

rule annotatePoNSV:
	input:
		svFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.txt",
		PoNFile="results/panelOfNormalsSV/PanelOfNormalsSV.txt",
		blackListFile="results/panelOfNormalsSV/PoNBlacklistBins.txt"
	output:
		outputSVAnnotFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe",
	params:
		annotScript=config["annotPoNSV_script"],
		svabafuncs=config["svaba_funcs"],
	log:
		"logs/combineSvabaGrocsvsTitan/{tumor}.annotPoNSV.log"
	shell:
		"Rscript {params.annotScript} --id {wildcards.tumor} --svaba_funcs {params.svabafuncs} --svFile {input.svFile} --PoNFile {input.PoNFile} --blackListFile {input.blackListFile} --outputSVAnnotFile {output.outputSVAnnotFile} 2> {log} > {log}"

rule filterSVs:
	input:
		svFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe"
	output:
		outputSVFiltFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.bedpe",
		outputSummaryFile="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.summary.txt"
	params:
		filterScript=config["filterSVs_script"],
		minFreqPoNSVBkptOverlap=config["PoN_minFreqSVbkpts"],
		# minFreqPoNCNVBkptOverlap=config["PoN_minFreqCNV"],
		minFreqPoNBlackList=config["PoN_minFreqBlackList"]
	log:
		"logs/combineSvabaGrocsvsTitan/{tumor}.filterSVs.log"
	shell:
		"Rscript {params.filterScript} --id {wildcards.tumor} --svFile {input.svFile} --minFreqPoNSVBkptOverlap {params.minFreqPoNSVBkptOverlap} --minFreqPoNBlackList {params.minFreqPoNBlackList} --outputSVFile {output.outputSVFiltFile} --outputSummary {output.outputSummaryFile} 2> {log} > {log}"

rule plotSvabaGrocsvsTitan:
	input:
		svabaVCF="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.bedpe",
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

rule plotCircos:
	input:
		svabaTitanBedpe="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.sv.annotPoN.bedpe",
		svabaTitanCN="results/combineSvabaGrocsvsTitan/{tumor}/{tumor}.svabaTitan.cn.txt"
	output:
		"results/plotCircos/{tumor}/{tumor}_Circos.pdf"
	params:
		plotCIRCOSscript=config["plotCircos_script"],
		genomeBuild=config["genomeBuild"]
	log:
		"logs/plotCircos/{tumor}/{tumor}_Circos.log"
	shell:
		"Rscript {params.plotCIRCOSscript} --id {wildcards.tumor} --svFile {input.svabaTitanBedpe} --cnFile {input.svabaTitanCN} --genomeBuild {params.genomeBuild} --outPlotFile {output} > {log} 2> {log}"


