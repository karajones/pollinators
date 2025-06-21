## A methodological comparison of eDNA derived from flowers and DNA derived from bulk samples of insects

> [!NOTE]
> This manuscript has been accepted for publication in *Molecular Ecology*. I will update with a link once it's been published. The scripts and data provided in this repository are intended to supplement the bioinformatics methods section.

## About the study

This study compared metabarcoding of environmental DNA (eDNA) from flowers with bulk mixed samples of arthopods taken from blue and yellow Japanese beetle vane traps using metabarcoding with three markers. Our goals were to: 
1. Compare richness estimates obtained from metabarcoding plant-derived eDNA versus metabarcoding DNA from bulk samples of arthropods (from vane traps) that were pulverized and homogenized
2. Investigate the impact of ASV clustering on arthropod detections and richness estimates
3. Evaluate the impact of primer bias on pollinator detections by comparing the taxa obtained from three metabarcoding primers (16S, COI, and Bombus)
4. Assess primer performance in the context of metabarcoding database completeness to determine how each factor may be affecting detections

## Description of scripts

These scripts were used to process sequences derived from the study. They are run roughly in order (though there may be a little back and forth needed between the BLAST and LULU curation steps). FastQ sequencing files that were used with these scripts are available as part of NCBI GenBank Project [PRJNA1189042](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1189042).

- `01_fastq_processing.sh` Executable (bash) shell script that renames fastq files, removes primers and small fragments using [Cutadapt](https://cutadapt.readthedocs.io/en/v3.5/index.html), and generates read statistics using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://github.com/MultiQC/MultiQC), and [Seqkit](https://bioinf.shenwei.me/seqkit/).
- `02_dada2_denoising.R` R script using [dada2](https://benjjneb.github.io/dada2/) to denoise trimmed reads and remove chimeras. The resulting merged amplicon sequence variants (ASVs) are filtered using [decontam](https://github.com/benjjneb/decontam) to remove contaminants found in blanks and size selected for amplicon length.
- `03_parse_BLAST_results.R` R script to filter and perform a lowest common ancestor (LCA) analysis on BLAST results.
- `04_ASV_curation_LULU.R` R script for curating ASVs using [LULU](https://github.com/tobiasgf/lulu) at four minimum match values and creation of [Phyloseq](https://joey711.github.io/phyloseq/) object with final data.
- `05_pollinator_analyses.R` Final processing and analysis of data.

<!--
## Description of data

`pollinators_bilsoda22.RData` is the final RData object created after running all of the above scripts. It contains the three Phyloseq objects (one for each marker: 16S, COI, and Bombus) combined into a single Phyloseq object named `ps_combined`. It includes an OTU table with read abundance, sample data table (i.e., metadata; see below for details), taxonomy table with taxonomic assignments, and DNAStringSet containing the sequences for each ASV. The ASV names have a prefix and underscore added to differentiate which marker they come from (e.g., COI_ASV1).

## Metadata
Brief description of sample variables included in the sample data class of the Phyloseq object (`pollinators_bilsoda22.RData`)
> [!IMPORTANT]
>  Not all variables were used in the study but are included here for completeness

| Variable | Description |
| --- | --- |
| ASV	| Unique identifier for amplicon sequence variant (ASV); prefix before underscore indicates marker used to detect the ASV |
| sample | Unique sample name |
| sample_type2	| Sample type: flowers (for eDNA from flowers), yellow or blue vane trap (for ground arthropods from vane traps) |
| date	| Date sample was collected |
| plot_ID	| Unique name for plot where sample was collected |
| temp	| Temperature in degrees Celsius on day sample was collected |
| wind	| Approximate wind speed on day sample was collected |
| sky	| Approximate cloud cover on day sample was collected |
| rain	| Whether it was raining on day sample was collected |
| flower	| Species of flower sample (or genus if flower could not be identified to species) |
| flower_confidence	| Level of confidence sample collector had in flower species identification |
| blank	| If TRUE, ASV is from a blank (i.e., no template control) |
-->

## Disclaimers

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government. Although this information product, for the most part, is in the public domain, it also may contain copyrighted materials as noted in the text. Permission to reproduce copyrighted items must be secured from the copyright owner.

Please see the license and disclaimer files for additional details.

