# Setting Up An EAGER Run

This page presumes you are roughly familiar what each step of a NGS pipeline does.

## Before you start

Make sure you know your how your libraries have been constructed (e.g. UDG treated, capture data?) and sequenced.

### File Organisation

All FASTQ files from the same sample should be in a sample specific directory. You should _not_ mix FASTQ files from different samples. You should also make sure all the file names are in the same format. This should follow the Illumina default of:

```
<SAMPLE_NAME>_S<X>_L00<X>_R<X>_00X.{fastq/fq}.gz
```
The files should share the same `<SAMPLE_NAME>` and `_S<X>_` (where X represents a number). The lane information `_L00<X>_` and read pairing `_R<X>_` can have differnet numbers, but must be in the same order. 

### Loading the EAGER GUI

If you are running EAGER from a server, and logging in via `ssh` remember to log into that server with `ssh -X` or `ssh -Y`. This is required to open the GUI windows.

### Modules

#### Input data

##### Select input *.fq/.fq.gz Files

This menu is aimed at describing the nature of your input data. You should make sure you know what your input files are - in terms of number of FASTQ files per sample, and what naming format they are in.

Changing options in this menu in _some_ cases switch in some default parameters downstream. 

* **Organism Type** Human is default. Switching to Bacteria or Other (i.e. animals) will turn on reporting of 'quality filtering' statistics in ReportTable.
* **Age of Dataset** To our knowledge, this doesn't make any changes.
* **Treated Data** To our knowledge this modifies settings within PMDtools, depending on whether you have performed no UDG treatment, half-UDG or (full-)UDG treatment on your libraries.
* **Pairment** This tells EAGER you have paired-end data. > If you're unsure if you have paired-end data, you can check this by looking at your raw FASTQ files and see if you have both `_R1_` _and_ `_R2_` files for each sample.
* **Capture data?** allows you to select a reference panel of SNPs from common aDNA SNP capture protocols that you may have performed on your libraries.
* **Calculate on target** This is used to tell EAGER to calculate the number reads and coverage of the positions listed in the capture data BED file. This is used to calculate SNP capture efficiency (when you compare with your previous shotgun data).
* **Input is already concatenated (skip merging)** To our knowledge, this is a semi-redundant option which _should_ turn off the clip-and-merging option. You should only use this when you have previously run your data through EAGER and want to _re-map_ previously clipped-and-merged data (i.e. not fresh from the sequencer, and has already been processed by AdapterRemoval or ClipandMerge)
* **Concatenate lanewise together** This option tells EAGER to combine, per sample, any files that are from the same read (`_R1_` or `_R2_`) but are from different lanes `_L001`, `_L002_`... (up to 8). If you're unsure whether you have multi-lane data, you can check this by looking at your raw FASTQ files and see if you have multiple files with different `_L00*_` names. 
  * You will _always_ have these files when you have NextSeq data. This options means the FASTQ files are literally stuck on to the ends of each other. 
  * When dealing with paired-end data, AdapterRemoval/ClipAndMerge looks for the same read IDs (i.e. cluster coordinate from the flow cell) to merge together so this is not an issue. 
  * However, you will **not** recieve FastQC sequencing quality assessments per lane. This makes it difficult to assess whether _different_ lanes had different sequencing efficenies.
* **MTCapture Data?** This is similar to capture data but for mitochondrial positions rather than SNP arrays (see above).

##### Select output folder

This is self explanatory. 

##### Select Reference

The input reference file must be uncompressed (not with `.gz` at the end) and end in `.fa` or `.fasta`. 

**Name of mitochondrial chromosome** gives you the option to get mapping statistics of a particular _entry_ in a multi-header `<REFERENCE>.fasta` file, in addition to the whole reference. 
  * You can see this if you have multiple lines starting with `>` in your `<REFERENCE>.fasta` file. 
  * This does not specifically have to be the mitochondrial chromosome nor a chromosome at all, although all information in the ReportTable output will call it 'MT'. 
  * If you do not have HG19 or GrCH37 human reference fastas, you need to put the first word of the header (i.e. until the first space and without the `>`) in the field. The GI numbers that window information refers to are deprecated IDs from NCBI and may not be in your `<REFERENCE>.fasta` file. 

> It is preferable, if you are going to run lots of EAGER runs in parallel (same time) rather than sequentially (one at a time), to pre-index your reference file before you set up the EAGER run. This can cause crashes if multiple EAGER runs try index the same file at the same time. You can do this by running the following three commands: `bwa index <REFERENCE>.fa`, `samtools faidx <REFERENCE>.fa` and `java -jar picard.jar CreateSequenceDictionary R=<REFERENCE>.fa O=<REFERNECE>.dict`.

#### Options

##### CPU Cores to be used

This is the maximum number of cores a module can use, if the module is multi-threaded. Multi-threaded means a particular calculation can be split up and run in parallel to speed up the process. The defaults are normally fine for this, unless you have particularly large sequencing data (e.g. a whole HiSeq lane), or very high endogenous DNA.

##### Memory in GB

This is the maximum amount of RAM (random-access-memory) a module can use. This is where the information during calculations are stored (but not 'written to disk' like files are. The defaults are normally fine for this, unless you have particularly large sequencing data (e.g. a whole HiSeq lane), or very high endogenous DNA.

##### Use system tmp dir

This option allows you to save temporary files in the computers `/tmp/` directory. If turned off, temporary directories are made within each EAGER module output directory. The default is fine here.

##### FastQC Analysis

FastQC is a tool that allows you assess the quality of the sequencing run. If you already know the sequencing run was sucessfull - i.e. if you are re-mapping your data - you can turn this off. However it is not a problem if left on as it uses very little hard-drive space (although will take longer for the EAGER run to finish).

##### Adapter RM / Merging

This module removes remaining adapters (i.e. the stretch of DNA that connects your DNA molecules to your indices), and if you have paired-end data, merges complementary reads from the `_R1_` and `_R2_` files. Merging increases the confidence of a particular base call if both bases in the forward and reverse reads are the same. It also removes low-quality bases from the end of reads, and removes reads if they go under a certain length.

There are two options: AdapterRemoval and ClipAndMerge. AdapterRemoval is more established and has since been shown to have better performance. 

Additional Options:
  * **Forward Read Adapter** The sequence of the forward adapter that the tool uses to compare against the read, to identify if the adapter is present. The default is a common Illumina sequence.
  * **Reverse Read Adapter** The sequence of the reverse adapter that the tool uses to compare against the read, to identify if the adapter is present. The default is a common Illumina sequence.
  * **Minimum Base Quality** The lowest base quality (i.e. the sequencer base call confidence) allowed at the end of a read to be retained. If the base quality goes lower, these bases will be removed from the read.
  * **Minimum Sequence Length** The minimum number of  bases a read needs to have to be retained. If the length goes under this, the read will be removed from the sequencing file.
  * **Minimum Adapter Overlap** The minimum number of bases the read adapter must overlap by, before the overlapped bases are removed from the read.
  * **Perform only adapter clipping** Reads are not merged and only adapter clipped. This is often only used for modern data where the insert size (i.e. DNA molecule is so long, the read pairs do not overlap.
  * **Keep merged only** This retains only paired-end reads which successfully overlapped and were merged. Any read that lost it's pair (e.g. the pair did not read the minimum sequence length), is discarded.

##### QualityFiltering

This performs a similar function to AdapterRemoval, but with reduced functionality. It is only needed for already AdapterRemoved and Collapsed data, e.g. from public data. You do not (actually, cannot) need to use this at the same time as AdapterRemoval. It allows you to remove reads shorter than a partcular length, and trim bases that have a lower than threshold base calling quality.

##### Mapping

Mapping the module is where your DNA reads are compared to your Reference genome and finds the best place the read aligns or matches. There are different algorithms here, however the the most common tool in ancient DNA is `bwa` which is designed historically (and serendipitously for aDNA) for short reads.

A few comments on the others: 

  * CircularMapper is a variant of BWA for circular genomes (such as mitochondrial or bacterial genomes), which allows reads to span across the start/end of a linear reference genome. I.e. it allows you to get even coverage across this 'join', where if you mapped purely linearly you would get a reduction in coverage because the reads either don't align enough or one or the other end of the reference is selected.
  * BWAMem is better for longer reads in modern DNA and has extra functionality.
  * BowTie2 is considered to be pretty equivalent to BWA and often comes down to personal preference
  * We have no experience of Stampy
  
###### BWA

* **Readgroup** a prefix given to the name read of each read given to the BAM file containing metadata about the sequencing run. You rarely need to change this.
* **BWA Seedlength (-l)** the length of each read used for 'seeding' when doing the fast look up of the read against the reference genome index before aligning the rest of the read. Human DNA people often turn this off (by setting e.g. 1024), but this makes run time very slow but with the benefit of getting more accurate alignment of reads.
* **BWA Max # diff (-n)** The number of mismatches (.e.g mutations) a seed can have when looking up the best place on the reference genome. A integer (i.e. without a decimal) is the exact number of mismatches a seed can have. A 'float' (i.e. with a decimal) scales the number of mismatches to the length of the read at different 'strengths'. This can be explored using the tool hosted [here](https://apeltzer.shinyapps.io/BWAmismatches/), written by [Alex Peltzer](https://github.com/apeltzer). The `bwa` default for modern data is 0.04. Ancient DNA researchers often relax this to account for damage and historical divergence from the reference genome. 
  * When dealing with non-UDG treated data, people of the Kircher-school set this to 0.01 (_see Shapiro, B., & Hofreiter, M. (Eds.). (2012). Ancient DNA: Methods and Protocols. Humana Press. https://doi.org/10.1007/978-1-61779-516-9_), whereas the Schubert school set this to 0.03 (_see Schubert, M., Ginolhac, A., Lindgreen, S., Thompson, J. F., Al-Rasheid, K. A. S., Willerslev, E., … Orlando, L. (2012). Improving ancient DNA read mapping against modern reference genomes. BMC Genomics, 13, 178. https://doi.org/10.1186/1471-2164-13-178_), however it doesn't make a huge amount of difference.
  * For full-UDG treated data this value is often made more strict and increased to 0.1.
* **BWA Qualityfilter (-q)** This will remove reads below a particular mapping quality score. 
  * If you wish to allow reads to map to two places on the genome equally well to be kept, set this to 0. This is useful when there are large repetitive regions in your genome such humans and other eukaryotes.
  * If you wish to have only reads that can map to a single place set this to '37' (the maximum value). This is often done in ancient pathogen and bacteria work.
* **Filter unmapped Reads** Keeping this on will remove all unmapped reads in your final BAM file. If it is turned off those reads will be retained.
* **Extracted Mapped/Unmapped Reads** If this is turned on, it will separate the unmapped reads from the mapped reads, and each will be stored in separate files.

###### CircularMapper

Being a variant of `bwa`, this algorithm has mostly the same parameters as above, but with the following differences:

Removed paramters:

* **Readgroup**
* **Filter unmapped Reads**

Additional parameters:

* **Elongation factor** the number of base pairs to be copied from the end of the reference genome, and appended to the beginning (and vice versa).
* **Reference to extend** which fasta entry in your reference file to apply the elongation into (without the `>`). 
 
###### BWAMem

There are no advanced parameters for this algorithm

###### Bowtie 2

You are able to supply additional parameters for this algorithm, however you must set these in the same way as you would set up the program on the command line. For this we recommend that you check the `bowtie2` [manual page](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

###### Stampy

There are no advanced parameters for this algorithm

#### Complexity Estimation

#### Remove Duplicates

#### PMDtools

#### Contamination Estimation

#### Damage Calculation

#### SNP Calling

#### SNP Filtering

#### VCF2Genome

#### Clean Up

#### Create Report?

#### Generate Config File

This button will then create per-input sample output directories, which inside of each a `.xml` config file will be generated, named with the date and time.

If you get a 'success' message - check your output directory you got the expected number sample directories and `.xml` files

### Extras

#### Setting up multiple runs in a single EAGER-GUI session

If you wish to set up multiple EAGER runs, for example the same input data to different reference files, you should note the following things to make setting up the next ones easier. Once the 'Generate Config File' module has successfully executed, the green font of the three input buttons turn orange. This orange font indicates that the information for those input settings have not been updated - i.e. not changed since the previous run set up. 

You can use this information to check you've successfully updated the corresponding input data you wish to change. Using the same example as above: 
 
 1. I set up the first run successfully
 2. The input button font colours change orange
 3. I select a new reference FASTA
 4. The font colour for the **Select Reference** button changes green
 5. Before I press **Generate Config File** I check the __Select input .fq/.fq.gz Files__ and **Select output folder** buttons are still orange (have not changed), but the **Select Reference** button is green
 6. If yes, I press **Generate Config Files**
 
 With this you can double check you've not accidently changed something you didn't plan to.
 
 **Important** In certain cases, the `VCF2Genome` module will for some reason turn itself on after generating the first run. Subsequent runs will then try running that module and may crash. Always double check both the module hasn't turned on in the GUI window, and that it is not turned on in the config `.xml` file. Check in the [debugging page](debugging.md) how to fix a config file
 

