# RNA-Seq Quality Control Adaptor Trimming

# ** IN PROGRESS !!!! **

## Introduction




## Recommended tool

There are several programs available to perform  adapter trimming, but one of the most flexible is [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). It's been specifically designed for use with Illumina data, and importantly supports combined trimming of both single-ended and paired-ended reads. Furthermore, the Trimmomatic downlaod contains FASTA files of the technology dependent adapter sequences from Illumina (which are usually only available from Illumina directly).


## How to use

You can download the latest version of Trimmomatic (0.33) from [here](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip). The Zip file contains the trimmomatic .jar file, along with a set of Illumina specific FASTA files `adapters` directory. To check everything;'s working, use `java -jar` to *run* `Trimmomatic` - it should print out some usage information:

	[ben@icmw368 ~]$ java -jar /opt/biosofts/Trimmomatic-0.32/trimmomatic-0.32.jar 
	  Usage: 
		PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   	  or: 
        	SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>...
	[ben@icmw368 ~]$ 

#### Side note: aliasing `trimmomatic`

To aid readability, you may find it useful to `alias` the `java -jar /PATH/TO/trimmomatic-0.33.jar` command to something less unwieldy. Running the following command sets this up for you:

	alias trimmomatic="java -jar /PATH/TO/trimmomatic-0.33.jar"

Then you can just type `trimmomatic <options>` whenever you need to run it.


### Single vs Paired end reads

You should specify SE (single ended) **or** PE (paired ended).  For single ended reads, you'll only have a single input file so your commandline will look like:

`trimmomatic   SE   input_fastq_name.fastq  output_trimmed_name.fastq   trim_module_1   trim_module_2 ....

All reads from `input_fastq_name.fastq` that "pass" the trimming modules (i.e. have sufficient quality, length, etc etc) will be written to `output_trimmed_name.fastq`. 

### Paired end trimming

However, things are a little more complicated for paired end (PE) invocations. In this case, you'll have two input files (forwards and reverse), but 4 possible output files. These are (in the order they're written on the command line):

   - forwared_paired
   - forward_unpaired
   - reverse_paired
   - reverse_unpaired
   
If a particular pair of reads (both forward and reverse) `pass` the trimming step(s), they'll be written to the `forward_paired` and `reverse_paired` files. If however only a single member of the pair survives the filtering, it will be written to the respective `forward_unpaired` or `reverse_unpaired file`. 


### Example command
An example will help to explain what's going on. Suppose that we have 3 pairs of reads to start with, read1, read2 and read3. The *forward* strand reads are stored in the `sample1_R1.fastq` file, and the *reverse* reads in `sample1_R2.fastq`. Also, let's further suppose that:

  - read1 is high quality in both forward and reverse.
  - read2 is high quality forward, but the reverse read suffered from a drop in quality early on (which is quite common with MiSeq paired end reads).
  - read3 is generally poor quality throughout.


We want to trim off any poor quality bases (say a base call quality is less than 20) at the ends of the reads, and only retain those reads that are at least 100 bases long after the quality trimming. We would do this with the following command line:

   `trimmomatic PE  sample1_R1.fastq sample1_R2.fastq  p_trim_sample1_R1.fastq u_trim_sample1_R1.fq p_trim_sample1_R2.fastq u_trim_sample1_R2.fastq    TRAILING:20   MINLEN:100`
   
The `TRAILING` module does the quality trimming (hence the `20`), and the `MINLEN` module ensures that any `post-trimmed` reads meat the minimum length requirements.

In terms of output, we would find that:

   - `p_trim_sample1_R1.fastq` contains the *forward* read for read1
   - `p_trim_sample1_R2.fastq` contains the *reverse* read for read1
   - `u_trim_sample1_R1.fastq` contains the *forward* read for read2
   - There is no entry for the *reverse* of read2 as it failed the trimming step
   - There are no entries at all for read3, as both *forward* and *reverse* were dropped by the filtering.



### Available modules

`trimmomatic` allows multiple trimming *operations*, which may be combined together or run one after the other. These range from simple "*cut n bases from the end of each read*" (TRAILING) to more complex methods to perform "*sliding window*" (SLIDINGWINDOW) and "*PCR primer/adapter*" (ILLUMINACLIP) trimming. 

The command for `trimmomatic` follows the general format:
   - `trimmomatic` - the command itself
   - `SE` or `PE` - work with single ended (`SE`) or paired-ended (`PE`)
   - `input_fastq_files

`trimmomatic   SE or PE   input_fastq_files   output_fastq_files  trimming_method_1 trimming_method_2 ....`

### Trimming methods

The general way of specifying each trimming method is NAME:argument1:argument2:... For the complete list of methods available in trimmomatic you should consult the [manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.33.pdf), but for getting started, you'll probably want to use some combination of the following:

   - **CROP** and **HEADCROP**
       - CROP and HEADCROP trim bases from the end and beginning of reads respectively, regardless of the base quality. In the case of CROP the single argument specifies the length of the resulting read (i.e. CROP:100 will filter reads such that they are at **most** 100 bases long), while HEADCROP's argument specifies the number of bases to remove from the *start* Thus HEADCROP:5 will remove 5 bases from the beginning of each read.
       
   - **MINLENGTH**
       - Drop any reads that are not at least `n` bases long. You would usually apply this method as the last step (after quality trimming, cropping etc). The single argument `n` specifies the mimum length of the read.
       
   - **SLIDINGWINDOW**
       - Trim the read once the average quality of a subset of bases (the *sliding window*) drops below a threshhold. The window *slides* from 5' to 3', thus trimming occurs at the *end* of the read.  Arguments are:
           - `windowSize` - how many base qualities should be considered in the averaging.
           - `requiredThreshold` - minimum average quality required before the read is trimmed.
           
   - **ILLUMINACLIP**
       - This (complex) method searches each read for any Illumina adapters present at either the start or end of the reads. If anything is found it's clipped off since it's biologically meaningless and can seriously confuse subsequent alignments and analysis. I would strongly recommend you read the documentation linked above to fully understand how ILLUMINACLIP works, and why it's important. It's particularly useful in the case of paired-end data. The arguments for ILLUMINACLIP are:
           - `adapter_file` - Fastq file containing the adapters to remove. There's a set of adapter files located in the  $$TRIMMOMATIC_ADAPTER_FILES directory on both the Linux and Mac machines. Otherwise you can specify your own.
           - `seed_mismatches` - How accurate must the initial adapter *match* be?
           - `plaindrome_clip_threshold` - See the documentation for the meaning of this.



