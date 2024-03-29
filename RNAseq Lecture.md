# RNAseq work flow using leture data

## Create the necessary working directory
```
cd $RNA_HOME
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
```
## Getting reference data (.fa file, in this case, chromosome 22)
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
## Copy the gene annotation files (.gtf file) to the working directory.
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

## Alignment using HISAT2
Perform alignments with HISAT2 to the genome and transcriptome.
#hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
 
 ```
hisat2 -p 8 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep1.sam
hisat2 -p 8 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep2.sam
hisat2 -p 8 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep3.sam

hisat2 -p 8 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep1.sam
hisat2 -p 8 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep2.sam
hisat2 -p 8 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep3.sam
```




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

## load the required libraries
```
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
```
## Load phenotype data from a file we saved in the current working directory
```
pheno_data = read.csv("UHR_vs_HBR.csv")
```
## Load ballgown data structure and save it to a variable "bg"
```
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
```
## Display a description of this object
```
bg
```
## Load all attributes including gene name
```
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
```
## Save the ballgown object to a file for later use
```
save(bg, file='bg.rda')
```

## Perform differential expression (DE) analysis with no filtering
```
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))
```

## Save a tab delimited file for both the transcript and gene results
```
write.table(results_transcripts, "UHR_vs_HBR_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

## Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
```
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
```

## Load all attributes including gene name
```
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])
```

## Perform DE analysis now using the filtered data
```
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))
```

## Output the filtered list of genes and transcripts and save to tab delimited files
```
write.table(results_transcripts, "UHR_vs_HBR_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

## Identify the significant genes with p-value < 0.05
```
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)
```

## Output the signifant gene results to a pair of tab delimited files
```
write.table(sig_transcripts, "UHR_vs_HBR_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "UHR_vs_HBR_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

## Exit the R session
```
quit(save="no")
```

## Back to bash
## What does the raw output from Ballgown look like?
```
head UHR_vs_HBR_gene_results.tsv
```

How many genes are there on this chromosome?
```
grep -v feature UHR_vs_HBR_gene_results.tsv | wc -l
```
How many passed filter in UHR or HBR?
```
grep -v feature UHR_vs_HBR_gene_results_filtered.tsv | wc -l
```
How many differentially expressed genes were found on this chromosome (p-value < 0.05)?
```
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | wc -l
```

## Display the top 20 DE genes
In the following commands we use grep -v feature to remove lines that contain “feature”. Then we use sort to sort the data in various ways. The k option specifies that we want to sort on a particular column (3 in this case which has the DE fold change values). The n option tells sort to sort numerically. The r option tells sort to reverse the sort

#Higher abundance in UHR
```
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t #Higher abundance in UHR
```

#Higher abundance in HBR by dropping r (reversing order)
```
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t #Higher abundance in HBR
```

Save all genes with P<0.05 to a new file.
```
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > DE_genes.txt
head DE_genes.txt
```
## ERCC DE Analysis
Compare the observed versus expected differential expression estimates for the ERCC spike-in RNAs:
```
cd $RNA_HOME/de/ballgown/ref_only
wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/assets/scripts/Tutorial_ERCC_DE.R
chmod +x Tutorial_ERCC_DE.R
./Tutorial_ERCC_DE.R $RNA_HOME/expression/htseq_counts/ERCC_Controls_Analysis.txt $RNA_HOME/de/ballgown/ref_only/UHR_vs_HBR_gene_results.tsv
```
results as pdf file


# edgeR Analysis
Make use of the raw counts you generated above using htseq-count
edgeR is a bioconductor package designed specifically for differential expression of count-based RNA-seq data
This is an alternative to using stringtie/ballgown to find differentially expressed genes
```
cd $RNA_HOME/
mkdir -p de/htseq_counts
cd de/htseq_counts
```
This next step creates a mapping file that will help us translate from ENSG IDs to Symbols. It does this by parsing the GTF transcriptome file we got from Ensembl.  That file contains both gene names and IDs. 
```
perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $RNA_REF_GTF | sort | uniq > ENSG_ID2Name.txt
head ENSG_ID2Name.txt
```
Determine the number of unique Ensembl Gene IDs and symbols. 
```
#count unique gene ids
cut -f 1 ENSG_ID2Name.txt | sort | uniq | wc -l

#count unique gene names
cut -f 2 ENSG_ID2Name.txt | sort | uniq | wc -l

#show the most repeated gene names
cut -f 2 ENSG_ID2Name.txt | sort | uniq -c | sort -r | head
```
Launch R:
```
R

```

###R code###
#Set working directory where output will go
```
working_dir = "~/workspace/rnaseq/de/htseq_counts"
setwd(working_dir)
```

#Read in gene mapping
```
mapping=read.table("~/workspace/rnaseq/de/htseq_counts/ENSG_ID2Name.txt", header=FALSE, stringsAsFactors=FALSE, row.names=1)
```

# Read in count matrix
```
rawdata=read.table("~/workspace/rnaseq/expression/htseq_counts/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
```

# Check dimensions
```
dim(rawdata)
```

# Require at least 25% of samples to have count > 25
```
quant <- apply(rawdata,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
rawdata <- rawdata[keep,]
dim(rawdata)
```


#################
# Running edgeR #
#################

# load edgeR
```
library('edgeR')
```

# make class labels rep will repeat the labet at designated times
```
class <- factor( c( rep("UHR",3), rep("HBR",3) ))
```

# Get common gene names
```
genes=rownames(rawdata)
gene_names=mapping[genes,1]
```

# Make DGEList object
```
y <- DGEList(counts=rawdata, genes=genes, group=class)
nrow(y)
```
# TMM Normalization
```
y <- calcNormFactors(y)
```

# Estimate dispersion
```
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
```
# Differential expression test
```
et <- exactTest(y)
```

# Print top genes
```
topTags(et)
```

# Print number of up/down significant genes at FDR = 0.05  significance level
```
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(y)[as.logical(de)]
```

# Output DE genes
# Matrix of significantly DE genes
```
mat <- cbind(
 genes,gene_names,
 sprintf('%0.3f',log10(et$table$PValue)),
 sprintf('%0.3f',et$table$logFC)
)[as.logical(de),]
colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")
```

# Order by log fold change
```
o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]
```
# Save table
```
write.table(mat, file="DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")
```

#To exit R type the following
```
quit(save="no")
```

Once we have run the edgeR tutorial, compare the sigDE genes to those saved earlier from cuffdiff:
This next step creates a mapping file that will help us translate from ENSG IDs to Symbols. 

```
cat $RNA_HOME/de/ballgown/ref_only/DE_genes.txt
cat $RNA_HOME/de/htseq_counts/DE_genes.txt
```

Pull out the gene IDs
```
cd $RNA_HOME/de/

cut -f 1 $RNA_HOME/de/ballgown/ref_only/DE_genes.txt | sort  > ballgown_DE_gene_symbols.txt
cut -f 2 $RNA_HOME/de/htseq_counts/DE_genes.txt | sort > htseq_counts_edgeR_DE_gene_symbols.txt
```

## Visualize overlap with a venn diagram
http://www.cmbi.ru.nl/cdd/biovenn/
http://bioinfogp.cnb.csic.es/tools/venny/

#To get the two gene lists you could use cat to print out each list in your terminal and then copy/paste.
```
cat ballgown_DE_gene_symbols.txt
cat htseq_counts_edgeR_DE_gene_symbols.txt
```

cd $RNA_HOME/de/ballgown/ref_only
R

#####R######
```
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

pheno_data = read.csv("UHR_vs_HBR.csv")
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

bg

bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

save(bg, file='bg.rda')

results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id")

write.table(results_transcripts, "UHR_vs_HBR_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

write.table(results_transcripts, "KO_vs_RE_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "KO_vs_RE_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

write.table(sig_transcripts, "UHR_vs_HBR_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "UHR_vs_HBR_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

grep -v feature UHR_vs_HBR_gene_results_sig.tsv | wc -l
```

##Display the top 20 DE genes
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t #Higher abundance in KO
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t #Higher abundance in RE

# Construct a heatmap showcasing the significantly DE genes
```
R
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ballgown)
```
#Set your output pdf name
```
pdf(file="UHR_vs_HBR.pdf")
```
#Set working directory where results files exist
```
working_dir = "~workspace/rnaseq/de/ballgown/ref_only"
setwd(working_dir)
```

#List the current contents of this directory
```
dir()
```

#Loading object: Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline 
```
load('bg.rda')
```
#View a summary of the ballgown object
```
bg
```

#Load gene names for lookup later in the tutorial
```
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
```
#Pull the gene_expression data frame from the ballgown object
```
gene_expression = as.data.frame(gexpr(bg))
```
#Set the columns for finding FPKM and create shorter names for figures
```
data_columns=c(1:6)
```
#Calculate the differential expression results including significance
```
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))
results_genes[,"de"] = log2(results_genes[,"fc"])
```
## Write a simple table of differentially expressed transcripts to an output file

#Each should be significant with a log2 fold-change >= 2
```
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]

sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]
```
## Create a heatmap to vizualize expression differences between the eight samples

#Define custom dist and hclust functions for use with heatmaps
```
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

main_title="sig DE Transcripts"
par(cex.main=0.8)
sig_genes_de=sig_tn_de[,"id"]
sig_gene_names_de=sig_tn_de[,"gene_name"]

data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,col=rev(heat.colors(75)))


dev.off()

quit()
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

To get the two gene lists you could use cat to print out each list in your terminal and then copy/paste.

```
cat ballgown_DE_gene_symbols.txt
cat htseq_counts_edgeR_DE_gene_symbols.txt
```
