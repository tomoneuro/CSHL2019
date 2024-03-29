Down load reference from UCSC genome browser
Go to Tool --> Table browser --> select genome tab, BED browser extensibnle data --> create output
Using TABIX 

# Getting necessary data for sequencing processing and analysis
## Check conda channel
```
conda config --get channels
conda config --add channels default
conda config --add channels ggd-genomics
conda config --add channels bioconda
conda config --add channels conda-forge
conda install ggd
```

## Get reference data from UCSC GenomeBrowser and Ensambl
```
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/*
```

Getting pfam
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ucscGenePfam.txt.gz
```
Using bedtools, will try to open downloaded files...
```
bedtools intersect -a Homo_sapiens.GRCh38.98.gtf.gz -b ucscGenePfam.txt.gz -wa
Error: unable to open file or unable to determine types for file ucscGenePfam.txt.gz

- Please ensure that your file is TAB delimited (e.g., cat -t FILE).
- Also ensure that your file has integer chromosome coordinates in the 
  expected columns (e.g., cols 2 and 3 for BED).
```
Comes back with errors, which are common.
 
## look at those errors and correct
```
zcat ucscGenePfam.txt.gz | head 10
```  
need to remove colum 1
```
zcat ucscGenePfam.txt.gz | cut -f 2- > ucscGenePfam.bed
```
```
bedtools intersect -a Homo_sapiens.GRCh38.98.gtf.gz -b ucscGenePfam.bed -wa
***** WARNING: File Homo_sapiens.GRCh38.98.gtf.gz has inconsistent naming convention for record:
1	havana	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";

***** WARNING: File Homo_sapiens.GRCh38.98.gtf.gz has inconsistent naming convention for record:
1	havana	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
```

# ANACONDA and Go Get Data (ggd)


## Using GDD to get the right data
```
ggd search will provide the receipe of how to get and install the right data
```
```
ggd search grch37 gtfmpfam
```
#then following "receipe" provided by ggd search will download, process (including making $1 to chr#, $2 start, $3 end, etc.) for standardization.
```
ggd install grch37-pfam-domains-ucsc-v1
```


