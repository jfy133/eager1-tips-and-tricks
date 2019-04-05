# Setting Up An EAGER Run

This page presumes you are roughly familiar what each step of a NGS pipeline does.

## Before you start

Make sure you know your how your libraries have been constructed (e.g. UDG treated, capture data?) and sequenced.

## File Organisation

All FASTQ files from the same sample should be in a sample specific directory. You should _not_ mix FASTQ files from different samples. You should also make sure all the file names are in the same format. This should follow the Illumina default of:

```
<SAMPLE_NAME>_S<X>_L00<X>_R<X>_00X.{fastq/fq}.gz
```
The files should share the same `<SAMPLE_NAME>` and `_S<X>_` (where X represents a number). The lane information `_L00<X>_` and read pairing `_R<X>_` can have differnet numbers, but must be in the same order. 

## Loading the EAGER GUI

If you are running EAGER from a server, and logging in via `ssh` remember to log into that server with `ssh -X` or `ssh -Y`. This is required to open the GUI windows.

## Input Files Window (Select input *.fq/.fq.gz Files)

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

## Output directory (Select output folder)

This is self explanatory. 

## Select Reference

The input reference file must be uncompressed (not with `.gz` at the end) and end in `.fa` or `.fasta`. 

**Name of mitochondrial chromosome** gives you the option to get mapping statistics of a particular _entry_ in a multi-header `<REFERENCE>.fasta` file, in addition to the whole reference. 
  * You can see this if you have multiple lines starting with `>` in your `<REFERENCE>.fasta` file. 
  * This does not specifically have to be the mitochondrial chromosome nor a chromosome at all, although all information in the ReportTable output will call it 'MT'. 
  * If you do not have HG19 or GrCH37 human reference fastas, you need to put the first word of the header (i.e. until the first space and without the `>`) in the field. The GI numbers that window information refers to are deprecated IDs from NCBI and may not be in your `<REFERENCE>.fasta` file. 

> It is preferable, if you are going to run lots of EAGER runs in parallel (same time) rather than sequentially (one at a time), to pre-index your reference file before you set up the EAGER run. This can cause crashes if multiple EAGER runs try index the same file at the same time. You can do this by running the following three commands: `bwa index <REFERENCE>.fa`, `samtools faidx <REFERENCE>.fa` and `java -jar picard.jar CreateSequenceDictionary R=<REFERENCE>.fa O=<REFERNECE>.dict`.

## CPU Cores to be used

This is the maximum number of cores a module can use, if the module is multi-threaded. Multi-threaded means a particular calculation can be split up and run in parallel to speed up the process. The defaults are normally fine for this, unless you have particularly large sequencing data (e.g. a whole HiSeq lane), or very high endogenous DNA.

## Memory in GB

This is the maximum amount of RAM (random-access-memory) a module can use. This is where the information during calculations are stored (but not 'written to disk' like files are. The defaults are normally fine for this, unless you have particularly large sequencing data (e.g. a whole HiSeq lane), or very high endogenous DNA.

## Use system tmp dir

This option allows you to save temporary files in the computers `/tmp/` directory. If turned off, temporary directories are made within each EAGER module output directory. The default is fine here.

## FastQC Analysis

FastQC is a tool that allows you assess the quality of the sequencing run. If you already know the sequencing run was sucessfull - i.e. if you are re-mapping your data - you can turn this off. However it is not a problem if left on as it uses very little hard-drive space (although will take longer for the EAGER run to finish).

## Adapter RM / Merging

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

## QualityFiltering

This performs a similar function to AdapterRemoval, but with reduced functionality. It is only needed for already AdapterRemoved and Collapsed data, e.g. from public data. You do not (actually, cannot) need to use this at the same time as AdapterRemoval. It allows you to remove reads shorter than a partcular length, and trim bases that have a lower than threshold base calling quality.

## Mapping

## Complexitiy Estimation

## Remove Duplicates

## PMDtools

## Contamination Estimation

## Damage Calculation

## SNP Calling

## SNP Filtering

## VCF2Genome

## Clean Up

## Create Report?
