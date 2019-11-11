# Single cell RNA seq

## Illumina sequencing output
FASTQC assessment of the seq quality

## Creating 10X-compatible FASTQ files with cellranger mkfastq (likely done by CORE)

```
cellranger mkfastq \
	--localcores=12 \
	--run=/path/to/basecalls/ \
	--samplesheet=/path/to/SampleSheet.csv \
  ```
  
  ## These are the files produced by Cellranger mkfastq from a NextSeq500 sequencing run:
  ```
  $ ls /path/to/fastqs/

SeqCourse2018-10XGEX-LPLard_S1_L001_I1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L003_I1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L001_R1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L003_R1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L001_R2_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L003_R2_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L002_I1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L004_I1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L002_R1_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L004_R1_001.fastq.gz
SeqCourse2018-10XGEX-LPLard_S1_L002_R2_001.fastq.gz  SeqCourse2018-10XGEX-LPLard_S1_L004_R2_001.fastq.gz
```
  

  ## Primary data analysis using cellranger count
  Here is how to run Cellranger in local mode, if you have only a single workstation computer or if you want to restrict the run to a single node on your cluster. I find this is slow but reliable:
  ```
  SAMPLE=SeqCourse2018-10XGEX-LPLard
TRANSCRIPTOME=/path/to/transcriptome/folder/


cellranger count \
  --id=$SAMPLE \
  --jobmode=local \
  --localcores=12 \
  --transcriptome=$TRANSCRIPTOME \
	--fastqs=/path/to/folder/containing/your/fastqs/ \
	--sample=$SAMPLE \
  
  ```
  
  ```
    	--fastqs=/path/to/fastq/folder1/,/path/to/fastq/folder2/
```

```
  $samtools view possorted_genome_bam.bam 19:5795690-5802672 | head -n 1
NB551387:106:HFL3VBGX9:3:23601:23610:10303	0	19	5795691	255	56M	*	0	
```

