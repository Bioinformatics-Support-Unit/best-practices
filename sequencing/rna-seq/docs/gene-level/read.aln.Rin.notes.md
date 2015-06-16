Introduction
============

After read quality control, the next step is the map these reads to a reference genome.

-splice junction aln VS non-splice awareness alignment

The choice of tools depends on the organnism being seqeunced i.e. if it is a prokaryotic species you don't have to worry about splice junctions, but on the other hand if it is from eukaryotes eg. human, then the aligner will need to be able to handle splice junctions.

In this document, we will focus on aligning RNAseq data from a human sample. There are number of RNAseq aligners avialable e.g. Bt, BWa, tophat, hisat etc. Here, we will describe how to use Star (ref). STAR 2-pass mapping procedure has been recommended as the most accurate tool for RNAseq [ref].

In brief, in the STAR 2-pass approach, splice junctions detected in a first alignment run are used to guide the final alignment.

procedure
---------
1) Download and install STAR (link)

STAR manual includes installation instruction.

2) Generate STAR genome index of the reference genome

Note that STAR provide a pre-compiled genome index for human and mouse genomes which can be downloaded from:

[STARgenomes Index](http://it-collab01.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/)


If you need to build a STAR genome from scratch, what you need to do is:

  - Obtain a reference genome in *fasta* format.
  - Build a STAR genome index using

  ```STAR --runMode genomeGenerate --genomeDir outputDir --genomeFastaFiles reference.fa```

  (see STAR manual for more information)

To run STAR 2-pass alignment steps:
  1. First round: align reads to the reference genome
  ```STAR --genomeDir genomeDir --readFilesIn read1.fastq.gz read2.fastq.gz --readFilesCommand zcat --outStd SAM```

Note that this step will output a SAM file, which can then be compressed into BAM format using ```samtools```


  2. Use the spliced junction info (file named with an extension *.SJ.out.tab*) file generated from Step 1 to create a new index

  ```STAR --runMode genomeGenerate --genomeDir genomeDir_SJ --genomeFastaFiles ref.fa --sjdbFileChrStartEnd SJ.out.tab --sjdbOverhang 75```

  3. Second round alignment: align reads to the reference index generated in Step 2

  ```STAR --genomeDir genomeDir_SJ --readFilesIn read1.fastq.gz read2.fastq.gz --readFilesCommand zcat --outStd SAM```


Checking alignment quality
--------------------------

The mapping quality of the data can be quickly checked by looking the final mapping statistics information provided by STAR in Log.final.out – 


```                                   Finished on |       Jun 15 21:14:53
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
`````
