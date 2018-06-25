# Pipeline to run RNA-seq data processing

This is based on the pipeline used by Myrsini Kaforou.

**Sample key**

First you need to make a key listing which fastq files are associated with each
sample. This can be pasted from the emails from Core for each lane of data 
received and looks like this:

```
FileName	              Sample Name
WTCHG_384296_201191	    S97
WTCHG_384296_202179	    S328
WTCHG_384296_203167	    S81
```

Name this file "sample_fastq_key.txt" and place it in your base directory.

The "1.MakeSampleKey.R" script takes this file and converts it into a key for
read alignment, i.e. lists all the fastq files for a particular sample.


**Configuration file**

To make your configuration file, copy the template configuration file and edit
appropriately. Specifically, you need to input your base directory, the location
of your fastq files, and your sample key. 

As this pipeline is updated, different mapping parameters will be controlled 
through this file.


**Running the pipeline**

To run the whole pipeline in one go, log in to rescomp and use:

sh RunPipeline.sh "/path/to/your/config/file/config.txt"

This will make all the necessary scripts within your specified base directory and
submit them to the cluster in turn.


**Mapping**

- Tool used: STAR aligner (https://github.com/alexdobin/STAR)
- Reference data required: Genome & Annotation
	- It is a splice-aware aligner so it requires the annotation file
-  Commands used:
	
**Making Indexes:**

```
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEFASTA --sjdbGTFfile GTF --sjdbOverhang 99

[--runMode]: “genomeGenerate” to generate the index
[--runThreadN]: multi-threading
[--genomeDir]: Directory to place the index
[--genomeFastaFiles]: FASTA file for the genome
[-sjdbGTFfile]: annotation file in GTF format
[-sjdbOverhang]: (direct quote from the manual) specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junction database. For 2x100bp, the value should be 99.
```

**Mapping:**
```
STAR --genomeDir $GENOME --runThreadN 6 --readFilesIn $FASTQF $FASTQR --outFileNamePrefix $PREFIX --outSAMtype BAM Unsorted --outFilterMismatchNoverLmax 0.04 --outSAMprimaryFlag OneBestScore --outReads Unmapped Fastx

[--genomeDir]: specify reference genome file
[--runThreadN]: multi-threading
[readFilesIn]: input files in FASTQ format
[--outFileNamePrefix]: prefix of the output file
[--outSAMtype]: output type (SAM or BAM format). “Unsorted” means the BAM file will have paired-end reads always adjacent.
[--outFilterMismatchNoverLmax]: speicify the number of mismatches. “0.04” means having 0.04*200 = 8 mismatches maximum allowed (for 2x100bp paired-end reads).
[--outSAMprimaryFlag]: “OneBestScore” means only one alignment with the best score is primary
[--outReads]: “Unmapped Fastx” means that the unmapped reads will be in FASTQ format.
```

**Counting**
- One tool: HTSeq-count (http://www-huber.embl.de/HTSeq/doc/overview.html)
- Reference data required: Annotation (GTF)
- Command used:

```
python htseq-count -t exon -r name -m union --stranded=yes -i gene_name -f bam $BAM $GFF > $COUNTS

[-t]: “exon” to count the reads overlapping exon
[-r]: “name” if the BAM file have been sorted by name. If “BAM Unsorted” option was used in STAR, this should be used.
[-m]: several options. “union” seems to be the consensus.
[--stranded]: stranded or unstranded RNA-seq data.
[-i]: “gene_name” will look for the gene_name within the annotation file
[-f]: “bam” specifying the input format.
```

- Another tool: featureCounts (http://bioinf.wehi.edu.au/featureCounts/)
- Reference data required: Annotation (GTF)
- Command used:

```
featureCounts -s 1 -p -t gene -g gene_name -a $GTF -o $COUNTS $BAM

[-s]: “1” if data is stranded, “0” if data is unstranded.
[-p]: specifying that the data is paired-end
[-t]: “exon” to count the reads overlapping exon
[-g]: “gene_name” will look for the gene_name within the annotation file
[-a]: specifying the annotation file (GTF)
```