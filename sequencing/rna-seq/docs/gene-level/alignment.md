# Introduction

After read quality control, the next step is the map these reads to a reference genome. In order to assess the number of reads mapping to a particular gene, it is necessary to find the location of each read in a reference genome.

There are many tools available for aligning high throughput sequencing data (DNASeq and RNASeq). Some tools are capable of aligning reads across splice junctions ("splice aware") for example [STAR](https://github.com/alexdobin/STAR), [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml), and [HISAT](http://ccb.jhu.edu/software/hisat/manual.shtml). Some aligners such as [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [BWA](http://bio-bwa.sourceforge.net/) are not designed to map reads that span slice junctions.

The choice of tool for aligning RNASeq data is therefore dependent upon the organism in question. For example if it is a prokaryotic transcriptome you don't have to worry about splice junctions, but on the other hand if the transcriptomic data is from eukaryotes e.g. human, then the aligner will need to be able to handle splice junctions.

<!---
Give examples: Not suitable for high throughput sequencing data - BLAST
-->

In this document, we will use [STAR](https://github.com/alexdobin/STAR) to align a paired-end RNAseq data. We will describe how to run *STAR 2-pass* mapping procedure which has been shown by [Pär G Engström et al.](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2722.html) as one of the most accurate tool for mapping transcript reads to a genome. In brief, the STAR 2-pass approach uses splice junctions detected in a first alignment run to guide the final alignment.

### Procedure

1. Download and install STAR (link)

  STAR manual includes installation instruction.

2. Generate STAR genome index of the reference genome

Note that STAR provide a pre-compiled genome index for human and mouse genomes which can be downloaded from:

[STARgenomes Index](http://it-collab01.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/)


If you need to build a STAR genome from scratch, what you need to do is:

  - Obtain a reference genome in *fasta* format.
  - Build a STAR genome index using

  ```
  STAR --runMode genomeGenerate --genomeDir outputDir --genomeFastaFiles reference.fa
  ```

(see STAR manual for more information)

To run STAR 2-pass alignment steps:

1. First round: align reads to the reference genome

  ```
  STAR --genomeDir genomeDir --readFilesIn mate_1.fastq.gz mate_2.fastq.gz
   --readFilesCommand zcat
  ```

  Note that this step will produce a SAM file *Aligned.out.sam*, which can then be compressed to BAM format using `samtools`

2. Use the spliced junction info (file named with an extension *.SJ.out.tab*) file generated from Step 1 to create a new index. Note that --sjdboverhang is the length of the "overhang" on each side of a splice junctions.

  ```
  STAR --runMode genomeGenerate --genomeDir genomeDir_SJ
  --genomeFastaFiles ref.fa --sjdbFileChrStartEnd SJ.out.tab --sjdbOverhang 75
  ```

3. Second round alignment: align reads to the reference index generated in Step 2

  ```
  STAR --genomeDir genomeDir_SJ --readFilesIn mate_1.fastq.gz mate_2.fastq.gz
   --readFilesCommand zcat --outStd SAM
  ```

### Checking alignment quality

The mapping quality of the data can be quickly checked by looking at the final mapping statistics information provided by STAR in the file *Log.final.out*. Example output is shown below:

```
                                    Finished on |       Jun 15 21:14:53
       Mapping speed, Million of reads per hour |       42.42

                          Number of input reads |       42340824
                      Average input read length |       200
                                    UNIQUE READS:
                   Uniquely mapped reads number |       34033289
                        Uniquely mapped reads % |       80.38%
                          Average mapped length |       197.71
                       Number of splices: Total |       16987784
            Number of splices: Annotated (sjdb) |       16914333
                       Number of splices: GT/AG |       16864158
                       Number of splices: GC/AG |       99631
                       Number of splices: AT/AC |       7369
               Number of splices: Non-canonical |       16626
                      Mismatch rate per base, % |       0.14%
                         Deletion rate per base |       0.01%
                        Deletion average length |       1.85
                        Insertion rate per base |       0.01%
                       Insertion average length |       1.44
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       6971849
             % of reads mapped to multiple loci |       16.47%
        Number of reads mapped to too many loci |       14994
             % of reads mapped to too many loci |       0.04%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       3.06%
                     % of reads unmapped: other |       0.06%
```
