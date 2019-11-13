1 Finding number of variant = counting line of vcf except header
```zless platinum-exome.vcf.gz | grep -v "^#" | wc -l
177675
```

2 How many (biallelic) SNPs are in the file? (Hint: awk's "length" option or bcftools)
```bcftools filter -i 'TYPE="snp"' platinum-exome.vcf.gz | grep -v "^#" | wc -l
160150

bcftools stats platinum-exome.vcf.gz | grep ^SN
162280

SN	0	number of samples:	4
SN	0	number of records:	177675
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	162280
SN	0	number of MNPs:	4220
SN	0	number of indels:	12782
SN	0	number of others:	727
SN	0	number of multiallelic sites:	3371
SN	0	number of multiallelic SNP sites:	297

(my initial answer)
bcftools view -m2 -M2 -v snps platinum-exome.vcf.gz | wc -l 
159941
```
3 How many transitions? (Hint: awk '$4=="C" && $5=="T"')
```
zcat platinum-exome.vcf.gz | awk '$4=="C" && $5=="T"' | wc -l
22955
zcat platinum-exome.vcf.gz | awk '$4=="T" && $5=="C"' | wc -l
 22781
zcat platinum-exome.vcf.gz | awk '$4=="A" && $5=="G"' | wc -l
23198

zcat platinum-exome.vcf.gz | awk '$4=="G" && $5=="A"' | wc -l
23336
```
 ... sam up every transitions
 Total 92270

4 How many transversions? (Hint: awk '$4=="A" && $5=="C"')
```
zcat platinum-exome.vcf.gz | awk '$4=="A" && $5=="C"' | wc -l
16130
zcat platinum-exome.vcf.gz | awk '$4=="C" && $5=="A"' | wc -l
5485

zcat platinum-exome.vcf.gz | awk '$4=="A" && $5=="T"' | wc -l
3674

zcat platinum-exome.vcf.gz | awk '$4=="T" && $5=="A"' | wc -l
3703

zcat platinum-exome.vcf.gz | awk '$4=="C" && $5=="G"' | wc -l
8283

zcat platinum-exome.vcf.gz | awk '$4=="G" && $5=="C"' | wc -l
8086
```
... sum up evry transversions

67583

5. Is the transition to transversion ratio what you expect?
1.36

6. How many high confidence (QUAL >30) variants are there (Hint: bcftools filter)?
```
bcftools filter -i 'QUAL>30' platinum-exome.vcf.gz | grep -v "^#" | wc -l
122046


#(my answer)
bcftools filter -i'%QUAL>30' platinum-exome.vcf.gz | wc -l
122134


https://samtools.github.io/bcftools/howtos/index.html
```
7. How many variants are present as a heterozygous genotype in 1 individual? (Hint: the AC attribute in VCF)
```
bcftools filter -i 'AC=1' platinum-exome.vcf.gz | grep -v "^#" | wc -l
40351

#(my answer)
zcat platinum-exome.vcf.gz | grep "AC=1" | wc -l
39306

zcat platinum-exome.vcf.gz | grep "AC=1" | grep -v "^#" | wc -l

Charlotte Reames 10:22 AM
zcat platinum-exome.vcf.gz | grep -E  "AC=*1*;AF=" | wc -l
38926
```
8. Which chromosome has the most SNPs?
```
bcftools query -i 'TYPE="snp"' -f'%CHROM\n' platinum-exome.vcf.gz | grep -v "^#" | sort | uniq -c 
  16030 1
   6524 10
   9313 11
   8063 12
   2920 13
   5410 14
   5537 15
   7444 16
   9564 17
   2443 18
  12872 19
  10546 2
   4101 20
   2281 21
   4568 22
   7792 3
   5797 4
   6497 5
   8283 6
   8032 7
   5528 8
   7058 9
   3448 X
     99 Y

bcftools query -i 'TYPE="snp"' -f'%CHROM\n' platinum-exome.vcf.gz | grep -v "^#" | sort | uniq -c | awk '{print $2"\t"$1}'
1	16030
10	6524
11	9313
12	8063
13	2920
14	5410
15	5537
16	7444
17	9564
18	2443
19	12872
2	10546
20	4101
21	2281
22	4568
3	7792
4	5797
5	6497
6	8283
7	8032
8	5528
9	7058
X	3448
Y	99

(my attempts)

bcftools query -f '%CHROM\n' platinum-exome.vcf.gz| wc -l

zless platinum-exome.vcf.gz | grep "^#=1" | wc -l
```

9. Which chromosome has the highest density of SNPs?



10. What do you guess for the sex of the four individuals? How would you go about this?
Using PEDDY software to determine based on Het/Hom+alt ratio (should be higher in female than male)
http://home.chpc.utah.edu/~u1138933/peddy/platinum_run3/CSHL_AdvSeqTech_2019.html

11. NA12878, NA12891, and NA12892 are a family. Who do you think are is the kid and who are the parents?
Using PEDDY software to determine based on ISB2 number suggestive for portion of variants shared between two individuals (should be higher in parent-child relationship and lower in unrelated (parents or unrelated) individuals.
http://home.chpc.utah.edu/~u1138933/peddy/platinum_run3/CSHL_AdvSeqTech_2019.html
12. How would you go about finding de novo mutations in the kid?
Using GEMINI

```
gemini qery -q "SELECT # FROM fakevariants" - header learnSQL2.db

gemini query -q "select chrom,start,end,ref,alt from variants" trio.trim.vep.denovo.db
```
Now see which ones are in coding regions
```
gemini query -q "select chrom,start,end,ref,alt,is_coding from variants where is_coding=1" trio.trim.vep.denovo.db
```
Now see which ones are loss of function?
```gemini query -q "select chrom,start,end,ref,alt,is_coding from variants where is_coding=1 and is_lof=1" trio.trim.vep.denovo.db```

"and" connects filtering criterion
```
gemini query -q "select chrom,start,end,ref,alt,is_coding from variants where is_coding=1" trio.trim.vep.denovo.db | wc -l
6229
gemini query -q "select chrom,start,end,ref,alt,is_coding,impact,gene from variants where is_coding=1 and is_lof=1" trio.trim.vep.denovo.db | wc -l
94
gemini query -q "select chrom,start,end,ref,alt,is_coding,impact,gene from variants where is_lof=1" trio.trim.vep.denovo.db | wc -l
146
```
We also need to check the quality, not taking the face value.
Quality associated with those variants phred scale
quality score of vcf files QUAL=-10-log(Pfalse) same as Phred score
Phred score>=30
```
gemini query -q "select chrom,start,end,ref,alt,is_coding,impact,gene from variants where is_lof=1 and qual>=30" trio.trim.vep.denovo.db | wc -l
141
```

Phred score>=100
```
gemini query -q "select chrom,start,end,ref,alt,is_coding,impact,gene from variants where is_lof=1 and qual>=100" trio.trim.vep.denovo.db | wc -l
138
```
Another consideration is the frequency (alternate allele frequency) of variants and diseases max_aaf_all in all population, ancestry and find maximum value
```
gemini query -q "select chrom,start,end,ref,alt,is_coding,impact,gene from variants where is_lof=1 and qual>=100 and max_aaf_all<0.0002941176471" trio.trim.vep.denovo.db
chr2	21236250	21236251	G	A	1	stop_gained	APOB
chr2	54001608	54001609	G	T	1	stop_gained	CHAC2
chr15	57953654	57953655	C	CA	1	frameshift_variant	GCOM1
chr17	39973373	39973374	C	T	1	stop_gained	FKBP10
chr22	30218065	30218066	T	C	0	splice_acceptor_variant	ASCC2

```
Two critical factors - frequency of variants and impact (loss of function). Low frequency is likely pathogenic.


