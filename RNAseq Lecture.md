## RNAseq work flow using leture data

## Create the necessary working directory
```
cd $RNA_HOME
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
```
## Getting data
```
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ls
```
## View the first 10 lines of this file. Why does it look like this?
```
head chr22_with_ERCC92.fa
```
## How many lines and characters are in this file? How long is this chromosome (in bases and Mbp)?
```
wc chr22_with_ERCC92.fa
```
## View 10 lines from approximately the middle of this file. What is the significance of the upper and lower case characters?
```
head -n 425000 chr22_with_ERCC92.fa | tail
```
## What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?
```
cat chr22_with_ERCC92.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
```
## View a list of all sequences in our reference genome fasta file.
```
grep ">" chr22_with_ERCC92.fa
```
## Annotation
## Copy the gene annotation files to the working directory.
```
echo $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
```


## Indexing
```
cd $RNA_REFS_DIR
hisat2_extract_splice_sites.py $RNA_REF_GTF > $RNA_REFS_DIR/splicesites.tsv
hisat2_extract_exons.py $RNA_REF_GTF > $RNA_REFS_DIR/exons.tsv
hisat2-build -p 8 --ss $RNA_REFS_DIR/splicesites.tsv --exon $RNA_REFS_DIR/exons.tsv $RNA_REF_FASTA $RNA_REF_INDEX
ls
```

## Obtain RNAseq data

```
echo $RNA_DATA_DIR
mkdir -p $RNA_DATA_DIR
cd $RNA_DATA_DIR
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls
```

## Enter the data directory and view the first two read records of a file (in fastq format each read corresponds to 4 lines of data)
```
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
```
## Identify the following components of each read: read name, read sequence, and quality string

How many reads are there in the first library? Decompress file on the fly with ‘zcat’, pipe into ‘grep’, search for the read name prefix and pipe into ‘wc’ to do a word count (‘-l’ gives lines)
```
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l
```
## Prealignment QC using FASTQC
```
cd $RNA_HOME/data
fastqc *.fastq.gz
```
## Run MultiQC on your fastqc reports to generate a single summary report across all samples/replicates.
```
cd $RNA_HOME/data
multiqc .
```
## Flexbar trim
## First, set up some directories for output
```
echo $RNA_DATA_TRIM_DIR
mkdir -p $RNA_DATA_TRIM_DIR
```
## Download necessary Illumina adapter sequence files.
```
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa
```

cd $RNA_HOME
```
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22


cd $RNA_DATA_TRIM_DIR
fastqc *.fastq.gz
```

## Alignment

## SAM to BAM Conversion Convert HISAT2 sam files to bam files and sort by aligned position
```
samtools sort -@ 8 -o UHR_Rep1.bam UHR_Rep1.sam
samtools sort -@ 8 -o UHR_Rep2.bam UHR_Rep2.sam
samtools sort -@ 8 -o UHR_Rep3.bam UHR_Rep3.sam
samtools sort -@ 8 -o HBR_Rep1.bam HBR_Rep1.sam
samtools sort -@ 8 -o HBR_Rep2.bam HBR_Rep2.sam
samtools sort -@ 8 -o HBR_Rep3.bam HBR_Rep3.sam
```
## Merge HISAT2 BAM files
```
cd $RNA_HOME/alignments/hisat2
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=UHR.bam INPUT=UHR_Rep1.bam INPUT=UHR_Rep2.bam INPUT=UHR_Rep3.bam
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=HBR.bam INPUT=HBR_Rep1.bam INPUT=HBR_Rep2.bam INPUT=HBR_Rep3.bam
```
## Count the alignment (BAM) files to make sure all were created successfully (you should have 8 total)
```
ls -l *.bam | wc -l
ls -l *.bam
```

## Before we can view our alignments in the IGV browser we need to index our BAM files. We will use samtools index for this purpose. For convenience later, index all bam files.
```
echo $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR
find *.bam -exec echo samtools index {} \; | sh
```


## Review in IGV


### Expression analysis with stringtie

```
cd $RNA_HOME/
mkdir -p expression/stringtie/ref_only/
cd expression/stringtie/ref_only/
```
## stringtie <aligned_reads.bam> [options]* Extraoptions specificed below: 
‘-p 8’ tells Stringtie to use eight CPUs
‘-G ' reference annotation to use for guiding the assembly process (GTF/GFF3)
‘-e’ only estimate the abundance of given reference transcripts (requires -G)
‘-B’ enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
‘-o’ output path/file name for the assembled transcripts GTF (default: stdout)
‘-A’ output path/file name for gene abundance estimates

```
stringtie -p 8 -G $RNA_REF_GTF -e -B -o HBR_Rep1/transcripts.gtf -A HBR_Rep1/gene_abundances.tsv $RNA_ALIGN_DIR/HBR_Rep1.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o HBR_Rep2/transcripts.gtf -A HBR_Rep2/gene_abundances.tsv $RNA_ALIGN_DIR/HBR_Rep2.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o HBR_Rep3/transcripts.gtf -A HBR_Rep3/gene_abundances.tsv $RNA_ALIGN_DIR/HBR_Rep3.bam


stringtie -p 8 -G $RNA_REF_GTF -e -B -o UHR_Rep1/transcripts.gtf -A UHR_Rep1/gene_abundances.tsv $RNA_ALIGN_DIR/UHR_Rep1.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o UHR_Rep2/transcripts.gtf -A UHR_Rep2/gene_abundances.tsv $RNA_ALIGN_DIR/UHR_Rep2.bam
stringtie -p 8 -G $RNA_REF_GTF -e -B -o UHR_Rep3/transcripts.gtf -A UHR_Rep3/gene_abundances.tsv $RNA_ALIGN_DIR/UHR_Rep3.bam
 ```
 
## What does the raw output from Stringtie look like? For details on the Stringtie output files refer to Stringtie manual (outputs section)
```
less -S UHR_Rep1/transcripts.gtf
```
## View transcript records only and improve formatting
```
grep -v "^#" UHR_Rep1/transcripts.gtf | grep -w "transcript" | column -t | less -S
```
## Limit the view to transcript records and their expression values (FPKM and TPM values)
```
awk '{if ($3=="transcript") print}' UHR_Rep1/transcripts.gtf | cut -f 1,4,9 | less
```

## Gene and transcript level expression values can also be viewed in these two files:
```
column -t UHR_Rep1/t_data.ctab | less -S

less -S -x20 UHR_Rep1/gene_abundances.tsv
```

## Create a tidy expression matrix files for the StringTie results. This will be done at both the gene and transcript level and also will take into account the various expression measures produced: coverage, FPKM, and TPM.

```
cd $RNA_HOME/expression/stringtie/ref_only/
wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/assets/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl
```

## Create a tidy expression matrix files for the StringTie results. First TPM as expression measure. For this we use perl script so we can put everything at once

./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='<input bam files> ' --transcript_matrix_file=transcript_tpm_all_samples.tsv (transcript output file) --gene_matrix_file=gene_tpm_all_samples.tsv (gene output file)

```
./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='HBR_Rep1,HBR_Rep2,HBR_Rep3,UHR_Rep1,UHR_Rep2,UHR_Rep3' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv
```
## next FPKM as expression measure
```
./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='HBR_Rep1,HBR_Rep2,HBR_Rep3,UHR_Rep1,UHR_Rep2,UHR_Rep3' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv
```
## then, coverage as expression measure
```
./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='HBR_Rep1,HBR_Rep2,HBR_Rep3,UHR_Rep1,UHR_Rep2,UHR_Rep3' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv
```

## Visualizing transcript in tpm
```
column -t transcript_tpm_all_samples.tsv | less -S
```

## Visualizing gene in tpm
```
column -t gene_tpm_all_samples.tsv | less -S
```
## HTSEQ-COUNT
Run htseq-count on alignments instead to produce raw counts instead of FPKM/TPM values for differential expression analysis

htseq-count basic usage:
htseq-count [options] <sam_file> <gff_file>

’–format’ specify the input file format one of BAM or SAM. Since we have BAM format files, select ‘bam’ for this option.
’–order’ provide the expected sort order of the input file. Previously we generated position sorted BAM files so use ‘pos’.
’–mode’ determines how to deal with reads that overlap more than one feature. We believe the ‘intersection-strict’ mode is best.
’–stranded’ specifies whether data is stranded or not. The TruSeq strand-specific RNA libraries suggest the ‘reverse’ option for this parameter.
’–minaqual’ will skip all reads with alignment quality lower than the given minimum value
’–type’ specifies the feature type (3rd column in GFF file) to be used. (default, suitable for RNA-Seq and Ensembl GTF files: exon)
’–idattr’ The feature ID used to identify the counts in the output table. The default, suitable for RNA-SEq and Ensembl GTF files, is gene_id.

```
cd $RNA_HOME/
mkdir -p expression/htseq_counts
cd expression/htseq_counts

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/UHR_Rep1.bam $RNA_REF_GTF > UHR_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/UHR_Rep2.bam $RNA_REF_GTF > UHR_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/UHR_Rep3.bam $RNA_REF_GTF > UHR_Rep3_gene.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/HBR_Rep1.bam $RNA_REF_GTF > HBR_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/HBR_Rep2.bam $RNA_REF_GTF > HBR_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/HBR_Rep3.bam $RNA_REF_GTF > HBR_Rep3_gene.tsv
```

## Merge results files into a single matrix for use in edgeR. The following joins the results for each replicate together, adds a header, reformats the result as a tab delimited file, and shows you the first 10 lines of the resulting file :

```
cd $RNA_HOME/expression/htseq_counts/
join UHR_Rep1_gene.tsv UHR_Rep2_gene.tsv | join - UHR_Rep3_gene.tsv | join - HBR_Rep1_gene.tsv | join - HBR_Rep2_gene.tsv | join - HBR_Rep3_gene.tsv > gene_read_counts_table_all.tsv
echo "GeneID UHR_Rep1 UHR_Rep2 UHR_Rep3 HBR_Rep1 HBR_Rep2 HBR_Rep3" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv | column -t
```

-grep -v "__" is being used to filter out the summary lines at the end of the files that ht-seq count gives to summarize reads that had no feature, were ambiguous, did not align at all, did not align due to poor alignment quality, or the alignment was not unique.

-awk -v OFS="\t" '$1=$1' is using awk to replace the single space characters that were in the concatenated version of our header.txt and gene_read_counts_table_all.tsv with a tab character. -v is used to reset the variable OFS, which stands for Output Field Separator. By default, this is a single space. By specifying OFS="\t", we are telling awk to replace the single space with a tab. The '$1=$1' tells awk to reevaluate the input using the new output variable.

## ERCC expression analysis

```
cd $RNA_HOME/expression/htseq_counts
wget http://genomedata.org/rnaseq-tutorial/ERCC_Controls_Analysis.txt
cat ERCC_Controls_Analysis.txt

wget https://github.com/griffithlab/rnabio.org/raw/master/assets/scripts/Tutorial_ERCC_expression.pl
chmod +x Tutorial_ERCC_expression.pl
./Tutorial_ERCC_expression.pl
cat $RNA_HOME/expression/htseq_counts/ercc_read_counts.tsv

wget https://github.com/griffithlab/rnabio.org/raw/master/assets/scripts/Tutorial_ERCC_expression.R
chmod +x Tutorial_ERCC_expression.R
./Tutorial_ERCC_expression.R ercc_read_counts

```
# Differential expression
## Ballgown DE Analyis
Use Ballgown to compare the tumor and normal conditions. 

```
mkdir -p $RNA_HOME/de/ballgown/ref_only/
cd $RNA_HOME/de/ballgown/ref_only/
```
## Perform UHR vs. HBR comparison, using all replicates, for known (reference only mode) transcripts:
First create a file that lists our 6 expression files, then view that file, then start an R session where we will examine these results:

```
printf "\"ids\",\"type\",\"path\"\n\"UHR_Rep1\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep1\"\n\"UHR_Rep2\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep2\"\n\"UHR_Rep3\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep3\"\n\"HBR_Rep1\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep1\"\n\"HBR_Rep2\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep2\"\n\"HBR_Rep3\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep3\"\n" > UHR_vs_HBR.csv
cat UHR_vs_HBR.csv

R
```

## Following will be performed in R

###R code###

# load the required libraries
```
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
```
# Load phenotype data from a file we saved in the current working directory
```
pheno_data = read.csv("UHR_vs_HBR.csv")
```
# Load ballgown data structure and save it to a variable "bg"
```
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
```
# Display a description of this object
```
bg
```
# Load all attributes including gene name
```
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
```
# Save the ballgown object to a file for later use
```
save(bg, file='bg.rda')
```

# Perform differential expression (DE) analysis with no filtering
```
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))
```

# Save a tab delimited file for both the transcript and gene results
```
write.table(results_transcripts, "UHR_vs_HBR_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
```
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
```

# Load all attributes including gene name
```
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])
```

# Perform DE analysis now using the filtered data
```
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))
```

# Output the filtered list of genes and transcripts and save to tab delimited files
```
write.table(results_transcripts, "UHR_vs_HBR_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

# Identify the significant genes with p-value < 0.05
```
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)
```

# Output the signifant gene results to a pair of tab delimited files
```
write.table(sig_transcripts, "UHR_vs_HBR_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "UHR_vs_HBR_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

# Exit the R session
```
quit(save="no")
```


nano ./bashrc

export RNA_HOME=~/workspace/rnaseq
export RNA_DATA_DIR=$RNA_HOME/data
export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
export RNA_REFS_DIR=$RNA_HOME/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REF_INDEX.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2
export _JAVA_OPTIONS=-Djavax.accessibility.assistive_technologies=
