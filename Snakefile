
##### Differential expression analysis for Lexogen 5' Quantseq RNAseq data #####
##### 2020-10-13 JCF #####

##### Sources: #####
# https://www.lexogen.com/quantseq-data-analysis/

from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.1.2")

configfile: "config.yaml"

def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()
    
rule target:
	input:
		expand("results/fastqc/{sample}_fastqc.zip", 
				sample=config["samples"])
    
rule fastqc:
	input:
		get_samples
	output:
		"results/fastqc/{sample}_fastqc.html",
		"results/fastqc/{sample}_fastqc.zip"
	threads: 4
	conda: "envs/fastqc.yaml"	
	shell:
		"fastqc --outdir results/fastqc --format fastq --threads {threads} {input}"
