
##### Differential expression analysis for Lexogen 5' Quantseq RNAseq data #####
##### 2020-10-13 JCF #####

##### Sources: #####
# https://www.lexogen.com/quantseq-data-analysis/

from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.1.2")

configfile: "config.yaml"

OUTDIR = config["prefix"]

def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()
    
rule target:
	""" Remember when using f string- wildcards must be double-bracketed!"""
	input:
		expand(f"{OUTDIR}/results/qualimap/{{sample}}/qualimapReport.html", 
				sample=config["samples"]),
		expand(f"{OUTDIR}/results/STAR/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
				sample=config["samples"]),
                expand(f"{OUTDIR}/results/multiqc/{{pre}}_multiqc_report.html",
                        	pre=config["prefix"])
#		f"{OUTDIR}/results/counts/all.tsv" 
    
rule fastqc:
	input:
		get_samples
	output:
		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.html",
		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.zip"
	threads: 4
	params: 
		outdir = f"{OUTDIR}/results/fastqc"
	log: f"{OUTDIR}/logs/fastqc/{{sample}}.log"
	conda: "envs/fastqc.yaml"	
	shell:
		"fastqc --outdir {params.outdir} --format fastq --threads {threads} {input}"

rule create_polyA:
	"""Need polyA fasta for trimming. PolyA  """
	output:
		"resources/polyA.fa.gz"
	threads: 1
	shell:
		"""
		printf ">poly_a\nAAAAAAAAAAAAAAAAAA" | gzip > {output}
		"""

rule bbduk_trim:
	""" Read orientation: Seq primer1, Read, (may read into poly-A tail if insert is short), Seq primer2 (not used) """
	""" sequencing primer, index """
	input:
		qc = f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.html",
		fq = get_samples,
		polyA = "resources/polyA.fa.gz"
	output:
		f"{OUTDIR}/results/fq_trim/{{sample}}_trim.fq.gz"
	params:	
		adapters = config["adaptor_path"]
	threads: 4
	log: f"{OUTDIR}/logs/bbduk/{{sample}}.log"
	conda: "envs/bbmap.yaml"
	shell:
		"""
		bbduk.sh in={input.fq} out={output} ref={input.polyA},{params.adapters} \
		k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 t={threads} 2> {log}
		"""

rule STAR_index:
	""" Create STAR genome index, also create a dummy file for input to next rule."""
	""" Marked dummy file as temp, so it will be deleted on finsihing pipeline run."""
	input:
		fa = config["genome"],
		gtf = config["gtf"]
	output:
		output = temp("STAR.ok")
	params: 
		overhang = config["read_length"] - 1
	threads: 8
	conda: "envs/STAR.yaml"
	log: "logs/STAR/genome_index.log"
	shell:
		"""
		STAR --runThreadN {threads} \
		--runMode genomeGenerate -genomeDir STAR_index --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} \
		--sjdbOverhang {params.overhang} 2> {log} &&
		touch {output}	
		"""

rule STAR_map:
	input:
		fq = f"{OUTDIR}/results/fq_trim/{{sample}}_trim.fq.gz",
		dummy = "STAR.ok"	
	output:
		f"{OUTDIR}/results/STAR/{{sample}}_Log.final.out",
		f"{OUTDIR}/results/STAR/{{sample}}_Aligned.sortedByCoord.out.bam"
	params:
		genome = config["genome"],
		out = f"{OUTDIR}/results/STAR/{{sample}}_"
	conda: "envs/STAR.yaml"
	log: f"{OUTDIR}/logs/STAR/{{sample}}.log"
	threads: 8
	shell:
		"""
		STAR --runThreadN {threads} --genomeDir GenomeDir --readFilesIn {input.fq} --readFilesCommand zcat \
		--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
		--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
		--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.out} \
		--quantMode GeneCounts 2> {log}
		"""

rule bam_index:
	input:
		f"{OUTDIR}/results/STAR/{{sample}}_Aligned.sortedByCoord.out.bam"
	output:
		f"{OUTDIR}/results/STAR/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
	conda: "envs/samtools.yaml"
	threads: 4
	shell:
		"samtools index {input} -@ {threads}"

rule qualimap:
	""" Params following: """
	""" https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/03_QC_STAR_and_Qualimap_run.html#qualimap"""
	input:
		f"{OUTDIR}/results/STAR/{{sample}}_Aligned.sortedByCoord.out.bam"
	output:
		f"{OUTDIR}/results/qualimap/{{sample}}/qualimapReport.html"
	params:
		gtf = config["gtf"],
		outdir =  f"{OUTDIR}/results/qualimap/{{sample}}"
	threads: 4
	conda: "envs/qualimap.yaml"
	shell:
		"""
		unset display
		qualimap rnaseq -outdir {params.outdir} \
		-a proportional -bam {input} \
		-p strand-specific-forward -gtf {params.gtf} \
		--java-mem-size=8G
		"""	

rule multiqc:
	input:
               [f"{OUTDIR}/results/qualimap/{sample}/qualimapReport.html" for sample in config["samples"]]
	output:
               expand(f"{OUTDIR}/results/multiqc/{{pre}}_multiqc_report.html", 
			pre=config["prefix"])
	params:
		indir = f"{OUTDIR}/results",
		name = config["prefix"],
		outdir =  f"{OUTDIR}/results/multiqc"
	conda: "envs/multiqc.yaml"
	shell:
               "multiqc {params.indir} --outdir {params.outdir} --title {params.name}"

#rule count_matrx:
#	input:
#		[f"{OUTDIR}/results/STAR/{sample}_ReadsPerGene.out.tab" for sample in config["samples"]] 
#	output:
#		f"{OUTDIR}/results/counts/all.tsv"
#	params:
#		samples=list(config["samples"])
#	conda: "envs/pandas.yaml"
#	script:
#		"scripts/count-matrix.py"


