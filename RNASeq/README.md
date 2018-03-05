# Summary of RNA-seq data from the Abu Dhabi study

### Pilot Experiment

Samples from 19 individuals were included. 10 of these were controls and 9 had
type 2 diabetes. 7 individuals had follow-up samples included and 12 did not.
For each patient and timepoint, 8 different samples were used:

- Paxgene

- PBMC

- CD8

- CD14

- CD14_PBS

- CD14_LPS

- CD14_MET

- CD14_MET_LPS

The files from Core are stored here: /well/jknight/kate/RNAseq/P160233-REX

The README file gives the identifiers of the samples (n=208 in total). Extended
sample information is in the [file PilotSampleInfoExpanded.txt](PilotData/PilotSampleInfoExpanded.txt). This does not
include some information available on the parents of these individuals.

Alignment was performed by Core Genomics using hisat. The outputs including bam 
files, an alignment summary and QC files for each sample are found in individual 
folders for each sample (SampleName_hisat). 

Deduplication and exclusion of multiple mapping reads was performed. Feature 
counts were then calculated using Subread using the file accepted_hits.dedup_nameSorted.bam
to give a file called featureCount_dedup.txt for each sample. Anntotation was
using /data/Genomes/GTF/human_glk_v37.ERCC.gtf

A matrix of unnormalised feature counts was made using the
[CountMatrix.R](PilotData/CountMatrix.R) script and saved under
/well/jknight/AbuDhabiRNA/Pilot-count-data.txt. 

Differential expression analysis was then performed using DESeq2
[Pilot_RNA_seq_results.Rmd](PilotData/Pilot_RNA_seq_results.Rmd)