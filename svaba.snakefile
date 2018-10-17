configfile: "config/config.yaml"
configfile: "config/samples.yaml"

import glob
def getLRFullPath(base, filename):
  return glob.glob(''.join([base, "/*/outs/", filename]))


rule svabaAll:
  input: 
  	#expand("results/svaba/{tumor}", tumor=config["pairings"]),
  	#expand("results/svaba/{tumor}/{tumor}.svaba.somatic.sv.vcf", tumor=config["pairings"]),
  	#expand("results/svaba/{tumor}/{tumor}.svaba.unfiltered.somatic.sv.vcf", tumor=config["pairings"]),
  	expand("results/barcodeRescue/{tumor}.bxOverlap.vcf", tumor=config["pairings"])

rule runSvaba:
	input:
		tum=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], config["bamFileName"]),
		norm=lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], config["bamFileName"])
	output:		
		outDir="results/svaba/{tumor}",
		#somaticFiltVCF="results/svaba/{tumor}/{tumor}.svaba.somatic.sv.vcf",
		somaticUnfiltVCF="results/svaba/{tumor}/{tumor}.svaba.unfiltered.somatic.sv.vcf",
		bps="results/svaba/{tumor}/{tumor}.bps.txt.gz"
	params:
		svabaExe=config["svaba_exe"],
		refGenome=config["refGenome"],
		numThreads=config["svaba_numThreads"],
		dbSNPindelVCF=config["svaba_dbSNPindelVCF"],
		mem=config["svaba_mem"],
		runtime=config["svaba_runtime"],
		pe=config["svaba_numCores"]
	log:
		"logs/svaba/{tumor}.log"
	shell:
		"{params.svabaExe} run -t {input.tum} -n {input.norm} -G {params.refGenome} -p {params.numThreads} -D {params.dbSNPindelVCF} -a {output.outDir}/{wildcards.tumor} > {log} 2> {log}"


rule barcodeRescue:
	input:
		tumBam=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], config["bamFileName"]),
		unfiltVCF="results/svaba/{tumor}/{tumor}.svaba.unfiltered.somatic.sv.vcf",
		bps="results/svaba/{tumor}/{tumor}.bps.txt.gz"
	output:
		"results/barcodeRescue/{tumor}.bxOverlap.vcf"
	params:
		bxRescueScript=config["bxRescue_script"],
		id="{tumor}",
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		chrs=config["chrs"],
		minMapQ=config["bxRescue_minMapQ"],
		minLength=config["bxRescue_minLength"],
		windowSize=config["bxRescue_windowSize"],
		minRead=config["bxRescue_minReadOverlapSupport"],
		mem=config["bxRescue_mem"],
		runtime=config["bxRescue_runtime"],
		pe=config["std_numCores"]		
	log:
		"logs/barcodeRescue/{tumor}.bxOverlap.log"
	shell:
		"Rscript {params.bxRescueScript} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --id {params.id} --tumBam {input.tumBam} --vcf {input.unfiltVCF} --bps {input.bps} --chrs \"{params.chrs}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --minMapQ {params.minMapQ} --minLength {params.minLength} --windowSize {params.windowSize} --minReadOverlapSupport {params.minRead} --outFile {output} > {log} 2> {log}"
		