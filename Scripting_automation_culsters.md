# Scripting, automation and academic clusters

Any computational work will take more than one line of code. This is where bash comes in.

# Bash scripting: simple but powerful
```
#!/bin/bash                                                                                                                                                                                                 
#Author: Eric Bogenschutz

reference="$(ggd get-files grch37-reference-genome-ensembl-v1 -s 'Homo_sapiens' -g 'GRCh37' -p 'grch37-reference-genome-ensembl-v1.fa')"

minimap2 -t 4 --MD -ax map-pb $reference $kid_fq | samtools view -hSb | samtools sort > bams/chr22_HG002.grch37.bam

samtools index bams/chr22_HG002.grch37.bam

sniffles --genotype -t 4 -m bams/chr22_HG002.grch37.bam -v vcfs/chr22_HG002_unsort.grch37.vcf
```
To run those in more samples, we need bash.
```
vim align_call
```
i to insert (edit)
esc, then :wq will save and close the vim editor

vim.rtorr.com - vim is perhaps the best program for terminal running commands

## Once you created bash, you can run multiple samples at once
```
#List samples in family
SAMPLES=('HG_002' 'HG_003' 'HG_004')
#Iterate over the three samples
for i in "${SAMPLES[@]}" 
do
reference="$(ggd get-files grch37-reference-genome-ensembl-v1 -s 'Homo_sapiens' -g 'GRCh37' -p 'grch37-reference-genome-ensembl-v1.fa')"

minimap2 -t 4 --MD -ax map-pb $reference $i | samtools view -hSb | samtools sort > bams/chr22_$i.grch37.bam

samtools index bams/chr22_$i.grch37.bam

sniffles --genotype -t 4 -m bams/chr22_$i.grch37.bam -v vcfs/chr22_$i_unsort.grch37.vcf
done
```

# Simple practice
```
#! /bin/bash
for b in ./*.bam;
do
     echo "$b"
done
~      
#! /bin/bash
for c in ./*.bai;
do
     echo "$c"
done
~      

~      
```
```
#! /bin/bash
v=0

if [ $v == 0 ]
then
    echo "yes"

else
    echo "no"
fi
~       
```
no
ubuntu@ip-172-31-4-4:~$ vim iftest

#! /bin/bash
v=1

if [ $v == 0 ]
then
    echo "yes"

else
    echo "no"
fi
~      
```


