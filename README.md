# awesome-bioinformatics

This is a collection of awesome bioinformatics resources.

# Table of Contents
<!-- TOC -->

- [awesome-bioinformatics](#awesome-bioinformatics)
- [Table of Contents](#table-of-contents)
- [RNA-Seq](#rna-seq)
    - [QC](#qc)
    - [Alignment](#alignment)
    - [Quantification](#quantification)
    - [differential expression analysis](#differential-expression-analysis)
    - [detecting splice junctions](#detecting-splice-junctions)
    - [exploring genomic data](#exploring-genomic-data)
- [ChIP-Seq](#chip-seq)

<!-- /TOC -->

# RNA-Seq

## QC
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - A quality control tool for high-throughput sequence data
- [MultiQC](https://multiqc.info/) - A tool for aggregating bioinformatics results across many samples
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) - A tool for removing adapter sequences from high-throughput sequencing reads
- [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) - A wrapper tool for Cutadapt and FastQC
- [FastP](https://github.com/OpenGene/fastp) - A tool for preprocessing paired-end reads
## Alignment
- [STAR](https://github.com/alexdobin/STAR) - Rapid Spliced Transcript Alignment to a Reference
- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) - A fast and sensitive alignment tool for RNA-seq data
- [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) - A tool for aligning sequencing reads to a reference genome
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - A tool for aligning sequencing reads to a reference genome
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) - A fusion gene finder for RNA-seq data
- [FusionCatcher](https://github.com/ndaniel/fusioncatcher) - A tool for detecting fusion genes in RNA-seq data
- [Fusion-Inspector](https://github.com/FusionInspector/FusionInspector) - A tool for detecting fusion genes in RNA-seq data
- [ChimeraSlayer](https://github.com/sfu-compbio/ChimeraSlayer) - A tool for detecting chimeric reads in RNA-seq data
- [Sailfish](https://github.com/kingsfordgroup/sailfish) - A tool for quantifying the expression of transcripts from RNA-seq data
- [Kallisto](https://pachterlab.github.io/kallisto/) - A tool for quantifying the expression of transcripts from RNA-seq data

## Quantification
- [HTSeq](https://htseq.readthedocs.io/en/master/) - A Python framework to work with high-throughput sequencing data
- [RSeQC](https://github.com/broadinstitute/rseqc) - A set of quality control tools for RNA-seq data
- [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) - A tool to count reads mapped to features in RNA-seq experiments
- [Salmon](https://combine-lab.github.io/salmon/) - A tool for quantifying the expression of transcripts from RNA-seq data
- [RSEM](https://deweylab.github.io/RSEM/) - A software package for estimating gene and isoform expression levels from RNA-seq data
- [Sailfish](https://github.com/kingsfordgroup/sailfish) - A tool for quantifying the expression of transcripts from RNA-seq data
- [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) - A tool for analyzing RNA-seq data
- [StringTie](https://ccb.jhu.edu/software/stringtie/) - A tool for transcript assembly and quantification from RNA-seq data
- [Rail-RNA](https://github.com/nellore/rail) - A Python package for analyzing RNA-seq data
- [Seqtk](https://github.com/lh3/seqtk) - A tool for processing sequences in the FASTA or FASTQ format
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - A tool for aligning sequencing reads to a reference genome
- [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) - A tool for aligning sequencing reads to a reference genome and 

## differential expression analysis
- [DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) - A package for differential expression analysis
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) - A package for differential expression analysis
- [limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html) - A package for differential expression analysis
- [voom](https://www.bioconductor.org/packages/release/bioc/html/voom.html) - A package for differential expression analysis
- [DESeq](https://www.bioconductor.org/packages/release/bioc/html/DESeq.html) - A package for differential expression analysis
- [DEXSeq](https://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html) - A package for differential expression analysis


## detecting splice junctions
- [BWA](http://bio-bwa.sourceforge.net/) - A tool for aligning sequencing reads to a reference genome
- [Samtools](http://www.htslib.org/) - A tool for manipulating and analyzing high-throughput sequencing data
- [GATK](https://software.broadinstitute.org/gatk/) - A toolkit for analyzing high-throughput sequencing data
- [BEDTools](https://bedtools.readthedocs.io/en/latest/) - A set of tools for working with genomic intervals
- [IGV](https://software.broadinstitute.org/software/igv/) - A visualization tool for interactive exploration of genomic data
- [UCSC Genome Browser](https://genome.ucsc.edu/) - A web-based platform for exploring genomic data
- [Ensembl](https://www.ensembl.org/) - A web-based platform for exploring genomic data
- [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) - A web-based platform for exploring genomic data
- [UCSC Genome Browser Track Hubs](https://genome.ucsc.edu/goldenPath/help/trackDb/trackDbHub.html) - A web-based platform for 
## exploring genomic data
- [UCSC Genome Browser Wiki](https://genome.ucsc.edu/wiki/Main_Page) - A web-based platform for exploring genomic data
- [ENCODE](https://www.encodeproject.org/) - A web-based platform for exploring genomic data
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) - A tool for sequence similarity searching
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - A tool for sequence similarity

# ChIP-Seq

- [ChIP-Seq analysis pipeline](https://github.com/nf-core/chipseq) - A Nextflow-based ChIP-Seq analysis pipeline
