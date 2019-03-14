# EAGER1 Output Intepretation

## Report Table

Your first point of call after an EAGER1 run has finished is the `Report_<name>.csv` or `Report_<name>.html`.

We will go through each column of this table, describing what each one means and what you _typically_ expect from ancient DNA.

**# of Raw Reads prior Clip & Merge (C&M)** This column is essentially the number of reads in your FASTQ file(s). It represents all the reads that could be demultiplexed with your particular index combination. The number of reads should be close to the number you requested for sequenced, or in the case of paired-end, double requested number. Data from: FASTQC

**# reads after C&M prior mapping** Data from: AdapterRemoval
**# of Merged Reads** Data from: AdapterRemoval
**# reads not attempted to map** Data from: AdapterRemoval
**# mapped reads prior RMDup** Data from: Samtools
**# of Duplicates removed** Data from: DeDup/MarkDuplicates
**Mapped Reads after RMDup** Data from: DeDup/MarkDuplicates
**Endogenous DNA (%)** Data from: ReportTable
**Endogenous DNA Cap (%)** Data from: ReportTable
**Cluster Factor** Data from: ReportTable
**Mean Coverage** Data from: QualiMap
**std. dev. Coverage** Data from: QualiMap
**Coverage >= 1X in %** Data from: QualiMap
**Coverage >= 2X in %** Data from: QualiMap
**Coverage >= 3X in %** Data from: QualiMap
**Coverage >= 4X in %** Data from: QualiMap
**Coverage >= 5X in %** Data from: QualiMap
**# of reads on mitochondrium**  Data from: QualiMap
**AVG Coverage on mitochondrium**  Data from: QualiMap
**MT/NUC Ratio**  Data from: MtNucRatio
**DMG 1st Base 3'**  Data from: MapDamage/Damageprofiler
**DMG 2nd Base 3'**  Data from: MapDamage/Damageprofiler
**DMG 1st Base 5'**  Data from: MapDamage/Damageprofiler
**DMG 2nd Base 5'**  Data from: MapDamage/Damageprofiler
**average fragment length**  Data from: MapDamage/Damageprofiler
**median fragment length**  Data from: MapDamage/Damageprofiler
**GC content in %** Data from Qualimap

## Files
