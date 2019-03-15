# EAGER1 Output Intepretation

## Report Table

Your first point of call after an EAGER1 run has finished is the `Report_<name>.csv` or `Report_<name>.html`.

We will go through each column of this table, describing what each one means and what you _typically_ expect from ancient DNA.

**# of Raw Reads prior Clip & Merge (C&M)** The number of reads in your FASTQ file(s). This column represents all the reads that could be demultiplexed with your particular index combination. The number of reads should be close to the number you requested for sequencing, or in the case of paired-end, double requested number. Data from: FASTQC.

**# reads after C&M prior mapping** The number of reads after merging. This number will be slightly lower than the # of Raw Reads, due to removal of short fragments (<30bp) and adapter dimers. For paired-end data, this number will be closer to the requested number of reads for sequencing, since some pairs of reads have now been merged into a single read. Data from: AdapterRemoval/ClipAndMerge.

**# of Merged Reads** The number of reads that have been successfully merged. For single-end data, this field will be `0`. Data from: AdapterRemoval/ClipAndMerge.

**# reads not attempted to map** The reads that were excluded prior to mapping. Most of these are reads smaller than 30bp, which cannot be confidetly mapped, and are therefore excluded. Also included are reads that have been filtered due to low base quality scores. Data from: AdapterRemoval/ClipAndMerge.

**# mapped reads prior RMDup** The number of reads that mapped to the reference genome, before the removal of PCR duplicates. Data from: Samtools.

**# of Duplicates removed** The number of PCR duplicates removed. A PCR duplicate is a read that matches another read with no mismatches. `MarkDuplicates` will compare reads that have the same starting position, and keep the longest read while discarding the rest. `DeDup` will take into account both starting and ending position of a read, and ony compare among perfectly identical reads. Data from: DeDup/MarkDuplicates.

**Mapped Reads after RMDup** The number of reads remaining after the removal of PCR duplicates. This number will be equal to the number of mapped reads minus the number of duplicates. Data from: DeDup/MarkrDuplicates.

**Endogenous DNA (%)** The percentage of reads after adapter removal that were successfully mapped to the reference genome. In other words, the ratio of `# mapped reads prior RMDup` to the `# reads after C&M prior mapping`. Data from: ReportTable.

**Endogenous DNA Cap (%)** When a bed file is provided, and the option `Calculate on Target` is selected, this column will provide the percentage of reads after adapter removal that were successfully mapped to the regions of the reference genome that overlap the regions in the bed file. Data from: ReportTable.

**Cluster Factor** The number of times a read has been sequenced ___on average___ in a library. This is calculated by taking the ratio of `# mapped reads prior RMDup` to `Mapped Reads after RMDup`, and can be used as an estimate of library complexity. Cluster factors below `2` indicate the library is not yet sequenced to exhaustion. Cluster factors above `2` indicate most of the unique reads in a library have already been sequenced, so further sequencing will come at diminishing returns. Often, beyond cluster factors of `6` further sequencing is ill-advised. It should be noted that a high cluster factor in screening phase does not necessarily mean bad preservation, as it could be caused by a low-quality library (e.g. an over-amplified library). Data from: ReportTable.

**Mean Coverage** Data from: QualiMap.
**std. dev. Coverage** Data from: QualiMap.
**Coverage >= 1X in %** Data from: QualiMap.
**Coverage >= 2X in %** Data from: QualiMap.
**Coverage >= 3X in %** Data from: QualiMap.
**Coverage >= 4X in %** Data from: QualiMap.
**Coverage >= 5X in %** Data from: QualiMap.
**# of reads on mitochondrium**  Data from: QualiMap.
**AVG Coverage on mitochondrium**  Data from: QualiMap.
**MT/NUC Ratio**  Data from: MtNucRatio.
**DMG 1st Base 3'**  Data from: MapDamage/Damageprofiler.
**DMG 2nd Base 3'**  Data from: MapDamage/Damageprofiler.
**DMG 1st Base 5'**  Data from: MapDamage/Damageprofiler.
**DMG 2nd Base 5'**  Data from: MapDamage/Damageprofiler.
**average fragment length**  Data from: MapDamage/Damageprofiler.
**median fragment length**  Data from: MapDamage/Damageprofiler.
**GC content in %** Data from Qualimap.

## Files
