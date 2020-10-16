# 3prime-tag-seq


Process reads from Lexogen QuantSeq FWD libraries, which use a 3' tag approach for RNAseq. This library prep results in reads that  are strand-specific on the forward strand. When sequenced single-end, the read starts with sequencing primer 1 and reads into the insert. Sequencing primer 2 is after the polyA tail. If insert size is small, the read one may read into the polyA tail. 

Currently following approach outlined on Lexogen's website: https://www.lexogen.com/quantseq-data-analysis/
along with some QC steps, which consists of:
1. FASTQC on raw reads for quality
2. Trimming of TruSeq adaptors and any potential polyA tails with bbduk from bbtools.
3. Mapping with STAR.
4. Assessing bam files with Qualimap and tying all QC together with MultiQC to output a html report. 

Conda is used to manage environments on a per rule basis, so be sure to deploy Snakemake with the --use-conda flag. 

![DAG](https://github.com/JamieCFreeman/3prime-tag-seq/blob/main/README_files/3prime-tag-seg.svg?raw=true)

