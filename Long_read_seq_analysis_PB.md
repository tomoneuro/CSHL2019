# Packbio long read seq analysis

## install packages and data using conda and ggd
```
conda install minimap2 sniffles samtools gsort
ggd install grch37-reference-genome-ensembl-v1
ggd install grch37-hg002-pb-chr22-giab-v1
ggd install grch37-hg003-pb-chr22-giab-v1
ggd install grch37-hg004-pb-chr22-giab-v1
```
## create working directories
```
mkdir bams
mkdir vcfs
```
## Let's process samples: get fastq file, get reference then put into minimap2 to align, then switch 
minimap2 -t 4 --MD -ax <option a and x (preset for map-pb)> map-pb $reference $kid_fq | samtools view -hSb <another options h+S+b> | samtools sort > bams/chr22_HG002.grch37.bam <output bam file>
samtools index bams/chr22_HG002.grch37.bam <sort the bam file>
sniffles --genotype -t 4 -m bams/chr22_HG002.grch37.bam -v vcfs/chr22_HG002_unsort.grch37.vcf <variant calling using sniffles>
  
1. hg002 
```
kid_fq="$(ggd get-files grch37-hg002-pb-chr22-giab-v1 -p "*fastq.gz")"
reference="$(ggd get-files grch37-reference-genome-ensembl-v1 -s 'Homo_sapiens' -g 'GRCh37' -p 'grch37-reference-genome-ensembl-v1.fa')"
minimap2 -t 4 --MD -ax map-pb $reference $kid_fq | samtools view -hSb | samtools sort > bams/chr22_HG002.grch37.bam
samtools index bams/chr22_HG002.grch37.bam
sniffles --genotype -t 4 -m bams/chr22_HG002.grch37.bam -v vcfs/chr22_HG002_unsort.grch37.vcf
```

2. hg003
```
dad_fq="$(ggd get-files grch37-hg003-pb-chr22-giab-v1 -p "*fastq.gz")"
reference="$(ggd get-files grch37-reference-genome-ensembl-v1 -s 'Homo_sapiens' -g 'GRCh37' -p 'grch37-reference-genome-ensembl-v1.fa')"
minimap2 -t 4 --MD -ax map-pb $reference $dad_fq | samtools view -hSb | samtools sort > bams/chr22_HG003.grch37.bam
samtools index bams/chr22_HG003.grch37.bam
sniffles --genotype -t 4 -m bams/chr22_HG003.grch37.bam -v vcfs/chr22_HG003_unsort.grch37.vcf
```

3. hg004
```
mom_fq="$(ggd get-files grch37-hg004-pb-chr22-giab-v1 -p "*fastq.gz")"
reference="$(ggd get-files grch37-reference-genome-ensembl-v1 -s 'Homo_sapiens' -g 'GRCh37' -p 'grch37-reference-genome-ensembl-v1.fa')"
minimap2 -t 4 --MD -ax map-pb $reference $mom_fq | samtools view -hSb | samtools sort > bams/chr22_HG004.grch37.bam
samtools index bams/chr22_HG004.grch37.bam
sniffles --genotype -t 4 -m bams/chr22_HG004.grch37.bam -v vcfs/chr22_HG004_unsort.grch37.vcf
```

## Check, if we got right vcf: by counting line number ignoring header
```
grep -v "^#" vcfs/chr22_HG002_unsort.grch37.vcf | wc -l
600

grep -v "^#" vcfs/chr22_HG002_unsort.grch37.vcf | head -1
14	19817186	0	TATACATATATATAGATATATATAGATATATATAGATATACATATATATAGATATATATAGATATATAGATATATAT	N	.	PASS	IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=14;END=19817247;ZMW=9;STD_quant_start=32.549962;STD_quant_stop=32.591410;Kurtosis_quant_start=0.524114;Kurtosis_quant_stop=-0.017782;SVTYPE=DEL;SUPTYPE=AL;SVLEN=-61;STRANDS=+-;RE=11;REF_strand=7,7;AF=0.44	GT:DR:DV	0/1:14:11

grep -v "^#" vcfs/chr22_HG002_unsort.grch37.vcf | head -2
14	19817186	0	TATACATATATATAGATATATATAGATATATATAGATATACATATATATAGATATATATAGATATATAGATATATAT	N	.	PASS	IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=14;END=19817247;ZMW=9;STD_quant_start=32.549962;STD_quant_stop=32.591410;Kurtosis_quant_start=0.524114;Kurtosis_quant_stop=-0.017782;SVTYPE=DEL;SUPTYPE=AL;SVLEN=-61;STRANDS=+-;RE=11;REF_strand=7,7;AF=0.44	GT:DR:DV	0/1:14:11
20	23968448	1	N	<DUP>	.	PASS	IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=20;END=23969233;ZMW=8;STD_quant_start=18.713631;STD_quant_stop=77.989743;Kurtosis_quant_start=-1.649924;Kurtosis_quant_stop=4.606668;SVTYPE=DUP;SUPTYPE=SR;SVLEN=785;STRANDS=-+;RE=10;REF_strand=2,4;AF=0.625	GT:DR:DV	0/1:6:10
```

## Moving forward to sort+bgzip, index, and merge
```
genome=https://raw.githubusercontent.com/jbelyeu/ggd-recipes/master/genomes/Homo_sapiens/GRCh37/GRCh37.genome
gsort vcfs/chr22_HG002_unsort.vcf $genome | bgzip -c > vcfs/chr22_HG002.grch37.vcf.gz
gsort vcfs/chr22_HG003_unsort.vcf $genome | bgzip -c > vcfs/chr22_HG003.grch37.vcf.gz
gsort vcfs/chr22_HG004_unsort.vcf $genome | bgzip -c > vcfs/chr22_HG004.grch37.vcf.gz

tabix vcfs/chr22_HG002.grch37.vcf.gz
tabix vcfs/chr22_HG003.grch37.vcf.gz
tabix vcfs/chr22_HG004.grch37.vcf.gz
```

## Now, let's see the results - visualization example
```
samplot plot -c 22 -s 18057764 -e 18060464 -b bams/chr22_HG002.grch37.bam -t DEL -o del_hg002_22_18057764_18060464.png
```
