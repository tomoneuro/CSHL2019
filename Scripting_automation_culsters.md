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
we can now test with imput 1
```
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

```
#!/bin/bash

for bam_file_name in ./*;
do
    if [ $bam_file_name == "./chr22_HG002.bam" ]
    then
        echo "$bam_file_name"
        echo -e "a few more things\n"
    else
        echo "Nope, can't do it"
    fi
done
```
For multiple files with name chr22_HG00 in common,
```
for bam_file_name in ./chr22_HG00*;
do
    if [ $bam_file_name == "./chr22_HG00*.bam" ]
    then
        echo "$bam_file_name"
        echo -e "a few more things\n"
    else
        echo "Nope, can't do it"
    fi
done
```
## The codes we created for one sample, if succeeds, can be put in bash to be run automatically, when multiple sapmles were put together. Get familiar with looping, conditional criteria, etc. Make argument - best practice
Set up multiple input like this
```

#! /bin/bash
number=$1

if [ $1 == 0 ]
then
    echo "yes"

else
    echo "no"
fi

then
bash iftest 4
no
bash iftest 0
yes
```
For more larger samples
```
$bash align_call.sh argument_1 argument_2 argument_3     

#!/bin/bash                                                                                                                                                                                                 
#Author: Eric Bogenschutz

reference=$1
$sample=$2
                                                                                                                                                                 

minimap2 -t 4 --MD -ax map-pb $reference $sample | samtools view -hSb | samtools sort > bams/$sample.grch37.bam

samtools index bams/$sample.grch37.bam

sniffles --genotype -t 4 -m bams/$sample.grch37.bam -v vcfs/$sample_unsort.grch37.vcf
```

# multiple samples to run - large server time, failure can be big loss and big monetary loss. RUN PILOT before running seriously.
