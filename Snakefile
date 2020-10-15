
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
				sample=config["samples"]),
		expand("results/bam/{sample}/Aligned.out.bam", 
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

#rule multiqc:
#	input:
#		[f"results/fastqc/{sample}_fastqc.zip" for sample in config["samples"]]
#	output:
#		"results/multiqc.html"
#	conda: "envs/multiqc.yaml"
#	shell:
#		"multiqc results/fastqc"

rule create_polyA:
	"""Need polyA fasta for trimming. PolyA  """
	output:
		"resources/polyA.fa.gz"
	threads: 1
	shell:
		"""
		printf ">poly_a\nAAAAAAAAAAAAAAAAAA" | gzip > {output}
		"""

rule trim:
	""" Read orientation: Seq primer1, Read, (may read into poly-A tail if insert is short), Seq primer2 (not used) """
	""" sequencing primer, index """
	input:
		qc = "results/fastqc/{sample}_fastqc.html",
		fq = get_samples,
		polyA = "resources/polyA.fa.gz"
	output:
		"results/fq_trim/{sample}_trim.fq.gz"
	params:	
		adapters = config["adaptor_path"]
	threads: 4
	conda: "envs/bbmap.yaml"
	shell:
		"""
		bbduk.sh in={input.fq} out={output} ref={input.polyA},{params.adapters} \
		k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 t={threads}
		"""

rule STAR_index:
	input:
		fa = config["genome"],
		gtf = config["gtf"]
	output:
		directory('Dmel_STAR')
	threads: 8
	conda: "envs/STAR.yaml"
	shell:
		"""
		mkdir {output} && 
		STAR --runThreadN {threads} \
		--runMode genomeGenerate -genomeDir {output} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} \
		--sjdbOverhang 100
		"""

rule map:
	input:
		fq = "results/fq_trim/{sample}_trim.fq.gz",
		index = "Dmel_STAR"
	output:
		"results/bam/{sample}/Aligned.out.bam"
	params:
		genome = config["genome"]
	conda: "envs/STAR.yaml"
	threads: 8
	shell:
		"""
		touch {output}
		"""
