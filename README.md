# 3prime-tag-seq


Process reads from Lexogen QuantSeq libraries, which use a 3' tag approach for RNAseq.

Currently following approach outlined on Lexogen's website: https://www.lexogen.com/quantseq-data-analysis/
along with some QC steps, which consists of:
1. FASTQC on raw reads for quality
2. Trimming of TruSeq adaptors and any potential polyA tails with bbduk from bbtools.
3. Mapping with STAR.
4. Assessing bam files with Qualimap and tying all QC together with MultiQC to output a html report. 

