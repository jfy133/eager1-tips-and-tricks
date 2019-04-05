# EAGER1 Output Intepretation

## Report Table

Your first point of call after an EAGER1 run has finished is the `Report_<name>.csv` or `Report_<name>.html`.

We will go through each column of this table, describing what each one means and what you _typically_ expect from ancient DNA.

**# of Raw Reads prior Clip & Merge (C&M)** The number of reads in your FASTQ file(s). This column represents all the reads that could be demultiplexed with your particular index combination. The number of reads should be close to the number you requested for sequencing, or in the case of paired-end, double requested number. Data from: FASTQC.

**# reads after C&M prior mapping** The number of reads after merging. This number will be slightly lower than the # of Raw Reads, due to removal of short fragments (<30bp) and adapter dimers. For paired-end data, this number will be closer to the requested number of reads for sequencing, since some pairs of reads have now been merged into a single read. Data from: AdapterRemoval/ClipAndMerge.

**# of Merged Reads** The number of reads that have been successfully merged. For single-end data, this field will be `0`. Data from: AdapterRemoval/ClipAndMerge.

**# reads not attempted to map** The reads that were excluded prior to mapping. Most of these are reads smaller than 30bp, which cannot be confidetly mapped and are therefore excluded. Also included are reads that have been filtered due to low base quality scores. Data from: AdapterRemoval/ClipAndMerge.

**# mapped reads prior RMDup** The number of reads that mapped to the reference genome, before the removal of PCR duplicates. Data from: Samtools.

**# of Duplicates removed** The number of PCR duplicates removed. A PCR duplicate is a read that matches another read with no mismatches. `MarkDuplicates` will compare reads that have the same starting position, and keep the longest read while discarding the rest. `DeDup` will take into account both starting and ending position of a read, and ony compare among perfectly identical reads. Data from: DeDup/MarkDuplicates.

**Mapped Reads after RMDup** The number of reads remaining after the removal of PCR duplicates. This number will be equal to the number of mapped reads minus the number of duplicates. Data from: DeDup/MarkrDuplicates.

**Endogenous DNA (%)** The percentage of reads after adapter removal that were successfully mapped to the reference genome. In other words, the ratio of `# mapped reads prior RMDup` to the `# reads after C&M prior mapping`. Data from: ReportTable.

**Endogenous DNA Cap (%)** When a bed file is provided, and the option `Calculate on Target` is selected, this column will provide the percentage of reads after adapter removal that were successfully mapped to the regions of the reference genome that overlap the regions in the bed file. Data from: ReportTable.

**Cluster Factor** The number of times a read has been sequenced ___on average___ in a library. This is calculated by taking the ratio of `# mapped reads prior RMDup` to `Mapped Reads after RMDup`, and can be used as an estimate of library complexity. Cluster factors below `2` indicate the library is not yet sequenced to exhaustion. Cluster factors above `2` indicate most of the unique reads in a library have already been sequenced, so further sequencing will come at diminishing returns. Often, beyond cluster factors of `6` further sequencing is ill-advised. It should be noted that a high cluster factor in screening phase does not necessarily mean bad preservation, as it could be caused by a low-quality library (e.g. an over-amplified library). Data from: ReportTable.

**Mean Coverage** The average fold coverage across the genome. I.e. the number of reads covering each position on the genome on average. Data from: QualiMap.

**std. dev. Coverage** The standard deviation in coverage across the genome. Small standard deviations in coverage means that reads are distributed more uniformly across the genome, while a high value points to uneven distribution of reads across the genome. Data from: QualiMap.

**Coverage >= 1X in %** What proportion of the genome is covered by at least 1 read. Data from: QualiMap.

**Coverage >= 2X in %** What proportion of the genome is covered by at least 2 reads. Data from: QualiMap.

**Coverage >= 3X in %** What proportion of the genome is covered by at least 3 reads. Data from: QualiMap.

**Coverage >= 4X in %** What proportion of the genome is covered by at least 4 reads. Data from: QualiMap.

**Coverage >= 5X in %** What proportion of the genome is covered by at least 5 reads. Data from: QualiMap.

**# of reads on mitochondrium** The number of reads that map to the fasta entry you selected in the "Select Reference" menu. Typically this refers to the mitochondrial chromosome, hence the column name. Data from: QualiMap.

**AVG Coverage on mitochondrium** The average fold coverage on the mitochondrial chromosome, or other fasta entry selected in the "Select Reference" menu. When this is the mitochondrial chromosome, this number is expected to be higher than the mean coverage on the autosomes, since there are many more copies of the mitochondrial DNA than the nuclear DNA. Data from: MtNucRatio.

**MT/NUC Ratio** The ratio of `AVG Coverage on mitochondrium` to `Mean Coverage`. Data from: MtNucRatio.

**DMG 1st Base 3'** The proportion of reads that contain a C->T transition in the first base on the 3' terminus of the read.This value is expected to be higher than "DMG 2nd Base 3'", and about equal to "DMG 1nd Base 5'" for double stranded libraries. It is recommended to look at the fragment misincorporation plots to **ensure** the presence of a smooth "smiley" for non-UDG-treated libraries, or the lack of damage after the third 3' and 5' bases in the case of UDG-half-treated libraries. Data from: MapDamage/Damageprofiler. 

**DMG 2nd Base 3'** The proportion of reads that contain a C->T transition in the second base on the 3' terminus of the read.
This value is expected to be lower than "DMG 1nd Base 3'", and about equal to "DMG 2nd Base 5'" for double stranded libraries. It is recommended to look at the fragment misincorporation plots to **ensure** the presence of a smooth "smiley" for non-UDG-treated libraries, or the lack of damage after the third 3' and 5' bases in the case of UDG-half-treated libraries. Data from: MapDamage/Damageprofiler. 

**DMG 1st Base 5'** The proportion of reads that contain a G->A transition in the first base on the 5' terminus of the read. This value is expected to be higher than "DMG 2nd Base 5'", and about equal to "DMG 1nd Base 3'" for double stranded libraries. It is recommended to look at the fragment misincorporation plots to **ensure** the presence of a smooth "smiley" for non-UDG-treated libraries, or the lack of damage after the third 3' and 5' bases in the case of UDG-half-treated libraries. Data from: MapDamage/Damageprofiler. 

**DMG 2nd Base 5'** The proportion of reads that contain a G->A transition in the second base on the 5' terminus of the read. This value is expected to be higher than "DMG 1nd Base 5'", and about equal to "DMG 2nd Base 3'" for double stranded libraries. It is recommended to look at the fragment misincorporation plots to **ensure** the presence of a smooth "smiley" for non-UDG-treated libraries, or the lack of damage after the third 3' and 5' bases in the case of UDG-half-treated libraries. Data from: MapDamage/Damageprofiler. 

**average fragment length**  The mean fragment length of mapped reads. This value is expected to be low for ancient DNA (between 30-75bp, depending on the age and preservation of the material). Data from: MapDamage/Damageprofiler.

**median fragment length**  The median fragment lenght of mapped reads. This value is expected to be low for ancient DNA (between 30-75bp, depending on the age and preservation of the material). Data from: MapDamage/DamageProfiler.

**GC content in %** The GC content across all mapped reads. This value should be close to the GC proportion of your reference genome. When working with NextSeq single-end sequenced data, values significantly higher than expected may stem from poly-G tails. This is because the machine reads the lack of a light as a guanine incorporation, therefore adding a poly-G tail to short fragments. If you believe you suffer from this issue after looking into your fastq files, you can utilise a tool like `fastp` to correct your raw data. Data from Qualimap.

## Files
### Introduction
### 0-FastQC
### 1-AdapClip
### 3-Mapper
### 4-Samtools
### 5-DeDup
`*.cleaned_rmdup.sorted.bam` vs `*.sorted.cleaned.bam`
### 6-QualiMap
### 7-DnaDamage
Within a subdirectory with the name of your `*_rmdup.sorted.bam` you will find the output of MapDamage/DamageProfiler. 
In that folder, a pdf fragment misincorporation plot can be found named `Fragmisincorporation_plot.pdf` if you used MapDamage, or `DamagePlot.pdf` if you used DamageProfiler.
You can also find a file containing the deamination damage misincorporations per terminal base for each end of the fragments in:
`3pGtoA_freq.txt` & `5pCtoT_freq.txt`. This data can be used to make your own Damage misincorporation plot, if you are so inclined.
