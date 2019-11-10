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

## Create a tidy expression matrix files for the StringTie results. First TPM as expression measure
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

column -t gene_tpm_all_samples.tsv | less -S
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
