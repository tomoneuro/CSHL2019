1 Finding number of variant = counting line of vcf except header
zless platinum-exome.vcf.gz | grep -v "^#" | wc -l

2 How many (biallelic) SNPs are in the file? (Hint: awk's "length" option or bcftools)
bcftools view -m2 -M2 -v snps platinum-exome.vcf.gz | wc -l 
159941

3 How many transitions? (Hint: awk '$4=="C" && $5=="T"')
 zcat platinum-exome.vcf.gz | awk '$4=="C" && $5=="T"' | wc -l
22955

4 How many transversions? (Hint: awk '$4=="A" && $5=="C"')
zcat platinum-exome.vcf.gz | awk '$4=="A" && $5=="C"' | wc -l
16130

5. Is the transition to transversion ratio what you expect?

6. How many high confidence (QUAL >30) variants are there (Hint: bcftools filter)?
bcftools filter -i'%QUAL>30' platinum-exome.vcf.gz | wc -l
122134


https://samtools.github.io/bcftools/howtos/index.html
How many variants are present as a heterozygous genotype in 1 individual? (Hint: the AC attribute in VCF)
Which chromosome has the most SNPs?
Which chromosome has the highest density of SNPs?

What do you guess for the sex of the four individuals? How would you go about this?
NA12878, NA12891, and NA12892 are a family. Who do you think are is the kid and who are the parents?
How would you go about finding de novo mutations in the kid?
