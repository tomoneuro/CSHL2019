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
```
# Alignment with HISAT2
```
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
```
stringtie -p 8 -G ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf -e -B -o KO_Rep1/transcripts.gtf -A KO_Rep1/gene_abundances.tsv ~/workspace/rnaseq/team_exercise/alignment/hisat/KO_Rep1.bam 
stringtie -p 8 -G ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf -e -B -o KO_Rep2/transcripts.gtf -A KO_Rep2/gene_abundances.tsv ~/workspace/rnaseq/team_exercise/alignment/hisat/KO_Rep2.bam 
stringtie -p 8 -G ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf -e -B -o KO_Rep3/transcripts.gtf -A KO_Rep3/gene_abundances.tsv ~/workspace/rnaseq/team_exercise/alignment/hisat/KO_Rep3.bam


stringtie -p 8 -G ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf -e -B -o RE_Rep1/transcripts.gtf -A RE_Rep1/gene_abundances.tsv ~/workspace/rnaseq/team_exercise/alignment/hisat/RE_Rep1.bam
stringtie -p 8 -G ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf -e -B -o RE_Rep2/transcripts.gtf -A RE_Rep2/gene_abundances.tsv ~/workspace/rnaseq/team_exercise/alignment/hisat/RE_Rep2.bam 
stringtie -p 8 -G ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf -e -B -o RE_Rep3/transcripts.gtf -A RE_Rep3/gene_abundances.tsv ~/workspace/rnaseq/team_exercise/alignment/hisat/RE_Rep3.bam
```
## the raw output from Stringtie
```
less -S KO_Rep1/transcripts.gtf
```

## View transcript records only and improve formatting
```
grep -v "^#" KO_Rep1/transcripts.gtf | grep -w "transcript" | column -t | less -S
```
## Limit the view to transcript records and their expression values (FPKM and TPM values)
```
awk '{if ($3=="transcript") print}' RE_Rep1/transcripts.gtf | cut -f 1,4,9 | less
```
## Create a tidy expression matrix files for the StringTie results. This will be done at both the gene and transcript level and also will take into account the various expression measures produced: coverage, FPKM, and TPM. 
```
cd ~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/

wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/assets/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl
```
## Create a tidy expression matrix files for the StringTie results. First TPM as expression measure. For this we use perl script so we can put everything at once
```
./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='KO_Rep1,KO_Rep2,KO_Rep3,RE_Rep1,RE_Rep2,RE_Rep3' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='KO_Rep1,KO_Rep2,KO_Rep3,RE_Rep1,RE_Rep2,RE_Rep3' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='KO_Rep1,KO_Rep2,KO_Rep3,RE_Rep1,RE_Rep2,RE_Rep3' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv
```
# 
```
cd ~/workspace/rnaseq/team_exercise
mkdir -p expression/htseq_counts
cd expression/htseq_counts
```

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ~/workspace/rnaseq/team_exercise/alignment/hisat/KO_Rep1.bam ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf > KO_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ~/workspace/rnaseq/team_exercise/alignment/hisat/KO_Rep2.bam ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf > KO_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ~/workspace/rnaseq/team_exercise/alignment/hisat/KO_Rep3.bam ~/workspace/rnaseq/team_exercise/references/chr12_Homo_sapiens.GRCh38.95.gtf > KO_Rep3_gene.tsv

# Gene expression level
## Q1. Based on your stringtie results, what are the top 5 genes with highest average expression levels across all knockout samples? What about in your rescue samples? How large is the overlap between the two sets of genes? (Hint: You can use R for this analysis)
Here’s some R code to start you off:
### load your data into R
```
exp_table=read.table('gene_tpm_all_samples.tsv', header=TRUE)

exp_table[,'mean_ko'] = apply(exp_table[,c(2:4)], 1, mean)
exp_table[,'mean_rescue'] = apply(exp_table[,c(5:7)], 1, mean)

exp_table[order(exp_table$mean_ko,decreasing=T)[1:5],]
exp_table[order(exp_table$mean_rescue,decreasing=T)[1:5],]
```
## Use ballgown to perform differential analysis followed by visualization of their results.
```
mkdir -p ~/workspace/rnaseq/team_exercise/de/ballgown/ref_only/
cd ~/workspace/rnaseq/team_exercise/de/ballgown/ref_only/

printf "\"ids\",\"type\",\"path\"\n\"KO_Rep1\",\"KO\",\"~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/KO_Rep1\"\n\"KO_Rep2\",\"KO\",\"~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/KO_Rep2\"\n\"KO_Rep3\",\"KO\",\"~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/KO_Rep3\"\n\"RE_Rep1\",\"RE\",\"~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/RE_Rep1\"\n\"RE_Rep2\",\"RE\",\"~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/RE_Rep2\"\n\"RE_Rep3\",\"RE\",\"~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/RE_Rep3\"\n" > KO_vs_RE.csv


cat KO_vs_RE.csv
"ids","type","path"
"KO_Rep1","KO","~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/KO_Rep1"
"KO_Rep2","KO","~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/KO_Rep2"
"KO_Rep3","KO","~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/KO_Rep3"
"RE_Rep1","RE","~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/RE_Rep1"
"RE_Rep2","RE","~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/RE_Rep2"
"RE_Rep3","RE","~/workspace/rnaseq/team_exercise/expression/stringtie/ref_only/RE_Rep3"


R
```
#####R######
```
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

pheno_data = read.csv("KO_vs_RE.csv")
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

bg
ballgown instance with 11734 transcripts and 6 samples
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

save(bg, file='bg.rda')

results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id")

write.table(results_transcripts, "KO_vs_RE_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "KO_vs_RE_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

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

write.table(sig_transcripts, "KO_vs_RE_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "KO_vs_RE_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

grep -v feature KO_vs_RE_gene_results_sig.tsv | wc -l
```

##Display the top 20 DE genes
grep -v feature KO_vs_RE_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t #Higher abundance in KO
grep -v feature KO_vs_RE_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t #Higher abundance in RE

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
pdf(file="KOvsRE_DE.pdf")
```
#Set working directory where results files exist
```
working_dir = "~/workspace/rnaseq/team_exercise/de/ballgown/ref_only"
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




export RNA_HOME=~/workspace/rnaseq
export RNA_DATA_DIR=$RNA_HOME/data
export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
export RNA_REFS_DIR=$RNA_HOME/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REF_INDEX.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2
export _JAVA_OPTIONS=-Djavax.accessibility.assistive_technologies=
