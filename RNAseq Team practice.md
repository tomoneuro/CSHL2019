# RNAseq Team practice 11/9-10/19

## Create the necessary working directory

```
cd $RNA_HOME
mkdir -p ~/workspace/rnaseq/team_exercise/data
cd ~/workspace/rnaseq/team_exercise/data
```

## Getting RNA data
```
mkdir -p ~/workspace/rnaseq/team_exercise/references
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_B/dataset.tar.gz
tar -xzvf dataset.tar.gz
```
## Getting adapter reference
```
mkdir -p ~/workspace/rnaseq/team_exercise/references
cd ~/workspace/rnaseq/team_exercise/references
```

## Getting reference for Adapter trimming
```
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa
```

## Reference fasta corresponding to your team's assigned chromosome (e.g. chr12)
```
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chr12.fa
```

## Obtain annotated reference gtf file corresponding to your team's assigned chromosome (e.g. chr6)
```
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chr6_Homo_sapiens.GRCh38.95.gtf
```

## Indexing references
```
cd /home/ubuntu/workspace/rnaseq/refs

hisat2_extract_splice_sites.py chr12_Homo_sapiens.GRCh38.95.gtf > splicesites.tsv
hisat2_extract_exons.py chr12_Homo_sapiens.GRCh38.95.gtf > exons.tsv
hisat2-build -p 8 --ss splicesites.tsv --exon exons.tsv chr12.fa ~/workspace/rnaseq/team_exercise/references/chr12_index
```

## Prealignment QC
```
cd $RNA_HOME/data

fastqc *.fastq.gz
```
## review your fastqc data, then multiqc
```
multiqc .
```
## Trimming ~/workspace/rnaseq/team_exercise/data/trimmed/
```
echo $RNA_DATA_TRIM_DIR

mkdir -p $RNA_DATA_TRIM_DIR
```
## Download illumina adapter sequence file
```
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/rnaseq/team_exercise/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads SRR10045016_1.fastq.gz  --reads2 SRR10045016_2.fastq.gz --target ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045016_1

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/rnaseq/team_exercise/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads SRR10045017_1.fastq.gz  --reads2 SRR10045017_2.fastq.gz --target ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045017_1

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/rnaseq/team_exercise/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads SRR10045018_1.fastq.gz  --reads2 SRR10045018_2.fastq.gz --target ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045018_1

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/rnaseq/team_exercise/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads SRR10045019_1.fastq.gz  --reads2 SRR10045019_2.fastq.gz --target ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045019_1

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/rnaseq/team_exercise/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads SRR10045020_1.fastq.gz  --reads2 SRR10045020_2.fastq.gz --target ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045020_1

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/rnaseq/team_exercise/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads SRR10045021_1.fastq.gz  --reads2 SRR10045021_2.fastq.gz --target ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045021_1
```

## Post trimming QC
```
fastqc *.fastq.gz

echo $RNA_ALIGN_DIR
mkdir -p $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR


hisat2 -p 8 --rg-id=KO_Rep1 --rg SM:KO --rg LB:KO_Rep1 --rg PL:ILLUMINA --rg PU:SRR10045016 -x ~/workspace/rnaseq/team_exercise/references/chr12_index --dta --rna-strandness RF -1 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045016_1_1.fastq.gz -2 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045016_1_2.fastq.gz -S ./KO_Rep1.sam

hisat2 -p 8 --rg-id=KO_Rep2 --rg SM:KO --rg LB:KO_Rep2 --rg PL:ILLUMINA --rg PU:SRR10045017 -x ~/workspace/rnaseq/team_exercise/references/chr12_index --dta --rna-strandness RF -1 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045017_1_1.fastq.gz -2 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045017_1_2.fastq.gz -S ./KO_Rep2.sam

hisat2 -p 8 --rg-id=KO_Rep3 --rg SM:KO --rg LB:KO_Rep3 --rg PL:ILLUMINA --rg PU:SRR10045018 -x ~/workspace/rnaseq/team_exercise/references/chr12_index --dta --rna-strandness RF -1 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045018_1_1.fastq.gz -2 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045018_1_2.fastq.gz -S ./KO_Rep3.sam

hisat2 -p 8 --rg-id=RE_Rep1 --rg SM:RE --rg LB:RE_Rep1 --rg PL:ILLUMINA --rg PU:SRR10045016 -x ~/workspace/rnaseq/team_exercise/references/chr12_index --dta --rna-strandness RF -1 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045019_1_1.fastq.gz -2 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045019_1_2.fastq.gz -S ./RE_Rep1.sam

hisat2 -p 8 --rg-id=RE_Rep2 --rg SM:KO --rg LB:RE_Rep2 --rg PL:ILLUMINA --rg PU:SRR10045017 -x ~/workspace/rnaseq/team_exercise/references/chr12_index --dta --rna-strandness RF -1 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045020_1_1.fastq.gz -2 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045020_1_2.fastq.gz -S ./RE_Rep2.sam

hisat2 -p 8 --rg-id=RE_Rep3 --rg SM:RE --rg LB:RE_Rep3 --rg PL:ILLUMINA --rg PU:SRR10045018 -x ~/workspace/rnaseq/team_exercise/references/chr12_index --dta --rna-strandness RF -1 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045021_1_1.fastq.gz -2 ~/workspace/rnaseq/team_exercise/data/trimmed/SRR10045021_1_2.fastq.gz -S ./RE_Rep3.sam
```
## Post alignment QC
```
fastqc *.sam

multiqc .
```
## path information
```
nano ~/.bashrc
```

## Sam to Bam
```
samtools sort -@ 8 -o KO_Rep1.bam KO_Rep1.sam
samtools sort -@ 8 -o KO_Rep2.bam KO_Rep2.sam
samtools sort -@ 8 -o KO_Rep3.bam KO_Rep3.sam
samtools sort -@ 8 -o RE_Rep1.bam RE_Rep1.sam
samtools sort -@ 8 -o RE_Rep2.bam RE_Rep2.sam
samtools sort -@ 8 -o RE_Rep3.bam RE_Rep3.sam
```
## Merge bam files
```
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=KO.bam INPUT=KO_Rep1.bam INPUT=KO_Rep2.bam INPUT=KO_Rep3.bam
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=RE.bam INPUT=RE_Rep1.bam INPUT=RE_Rep2.bam INPUT=RE_Rep3.bam
```
## Indexing bam files
```
find *.bam -exec echo samtools index {} \; | sh
```

## Visalizing alignment and expression in IGV

## Use Stringtie to generate expression estimates from the SAM/BAM files 




###########DNA way####################
## Query name sort bam files
# Runtime: ~ 4min
java -Xmx60g -jar $PICARD SortSam I=KO.bam O=KO_namesorted.bam SO=queryname
java -Xmx60g -jar $PICARD SortSam I=RE.bam O=RE_namesorted.bam SO=queryname


## Mark duplicates
java -Xmx60g -jar $PICARD MarkDuplicates I=KO_namesorted.bam O=KO_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=KO_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
java -Xmx60g -jar $PICARD MarkDuplicates I=RE_namesorted.bam O=RE_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=RE_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT



## Position sort bam file
java -Xmx60g -jar $PICARD SortSam I=KO_namesorted_mrkdup.bam O=KO_sorted_mrkdup.bam SO=coordinate
java -Xmx60g -jar $PICARD SortSam I=RE_namesorted_mrkdup.bam O=RE_sorted_mrkdup.bam SO=coordinate

## Create bam index for use with GATK, IGV, etc

java -Xmx60g -jar $PICARD BuildBamIndex I=KO_sorted_mrkdup.bam
java -Xmx60g -jar $PICARD BuildBamIndex I=RE_sorted_mrkdup.bam
###########DNA way####################


export RNA_HOME=~/workspace/rnaseq
export RNA_DATA_DIR=$RNA_HOME/data
export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
export RNA_REFS_DIR=$RNA_HOME/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REF_INDEX.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2
export _JAVA_OPTIONS=-Djavax.accessibility.assistive_technologies=
