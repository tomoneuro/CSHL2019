# Using bedtools

First you need to sort files in order of chromosome#
```
 sort -k1,1 -k2,2n in.bed > in.sorted.bed
```
## Puzzles to help teach you more bedtools.
1. Create a BED file representing all of the intervals in the genome that are NOT exonic.
```
bedtools complement -i exons.bed -g genome.txt > non-exonic.bed
```
2. What is the average distance from GWAS SNPs to the closest exon? (Hint - have a look at the closest tool.)
```
bedtools closest -io -d -a gwas.bed -b exon.bed | head
bedtools closest -d -a gwas.bed -b exons.bed | awk '{ total += $11; count++ } END { print total/count }'
46713.1
```

alternative
```
bedtools closest -d -a exons.bed -b gwas.bed | awk '{ total += $11; count++ } END { print total/count }'
alternative
bedtools closest -a exons.bed -b gwas.bed -d | awk '{sum+=$11} END { print "Average = ",sum/NR}'
```

3. Count how many exons occur in each 500kb interval (“window”) in the human genome. (Hint - have a look at the makewindows tool.)
```
bedtools makewindows -g genome.txt -w 500000 > 500kwindow.txt
bedtools intersect -a 500kwindow.txt -b exons.bed | wc -l

ubuntu@ip-172-31-4-4:~/workspace/monday/bedtools$ bedtools makewindows -g genome.txt -w 500000 | head
chr1	0	500000
chr1	500000	1000000
chr1	1000000	1500000
chr1	1500000	2000000
chr1	2000000	2500000
chr1	2500000	3000000
chr1	3000000	3500000
chr1	3500000	4000000
chr1	4000000	4500000
chr1	4500000	5000000
ubuntu@ip-172-31-4-4:~/workspace/monday/bedtools$ bedtools makewindows -g genome.txt -w 500000 > 500kwindow.txt
ubuntu@ip-172-31-4-4:~/workspace/monday/bedtools$ bedtools intersect -a 500kwindow.txt -b exons.bed | head
chr1	11873	12227
chr1	12612	12721
chr1	13220	14409
chr1	14361	14829
chr1	14969	15038
chr1	15795	15947
chr1	16606	16765
chr1	16857	17055
chr1	17232	17368
chr1	17605	17742
ubuntu@ip-172-31-4-4:~/workspace/monday/bedtools$ bedtools intersect -a 500kwindow.txt -b exons.bed | wc -l
459987
ubuntu@ip-172-31-4-4:~/workspace/monday/bedtools$ bedtools makewindows -g genome.txt -w 500000 | wc -l
6343

45998/6343=72.5
```


4. Are there any exons that are completely overlapped by an enhancer? If so, how many?
```
grep -i "enhancer" hesc.chromHmm.bed > enhancer.bed
bedtools intersect -a exons.bed -b enhancer.bed -f 1 | head
bedtools intersect -a exons.bed -b enhancer.bed -f 1 | wc -l
13746
```

5. What fraction of the GWAS SNPs are exonic?
```
bedtools intersect -a gwas.bed -b non-exonic.bed | wc -l
16055
bedtools intersect -a gwas.bed -b exonic.bed | wc -l
3439
exonic 0.176
```

6. What fraction of the GWAS SNPs are lie in either enhancers or promoters in the hESC data we have?
```
grep -i "promoter" hesc.chromHmm.bed > promoter.bed
bedtools intersect -a gwas.bed -b promoter.bed | wc -l
404
bedtools intersect -a gwas.bed -b enhancer.bed | wc -l
881
cat gwas.bed | wc -l
17680

(404+881)/17680 = 0.0727
```

7. Create intervals representing the canonical 2bp splice sites on either side of each exon (don’t worry about excluding splice sites at the first or last exon). (Hint - have a look at the flank tool.)
```
bedtools flank -i exons.bed -g genome.txt -b 2| head
```
8. What is the Jaccard statistic between CpG and hESC enhancers? Compare that to the Jaccard statistic between CpG and hESC promoters. Does the result make sense? (Hint - you will need grep).
```
bedtools jaccard -a cpg.bed -b enhancer.bed
1148180	132977386	0.0086344	4969

bedtools jaccard -a cpg.bed -b promoter.bed
15661111	53551816	0.292448	20402
```
9. What would you expect the Jaccard statistic to look like if promoters were randomly distributed throughout the genome? (Hint - you will need the shuffle tool.)
bedtools shuffle -i <(grep Promoter hesc.chromHmm.bed) -g genome.txt \
  | sort -k1,1 -k2,2n \
> promoters.shuffled.bed

bedtools jaccard -a cpg.bed -b promoters.shuffled.bed
intersection	union	jaccard	n_intersections
341707	68521738	0.00498684	842

(wrong)
shuf promoter.bed > shuffledpromoter.bed
sort -k 1,1n shuffledpromoter.bed > shuffledpromoter1.bed

bedtools jaccard -a cpg.bed -b shuffledpromoter1.bed

bedtools shuffle -chrom -i promoter.bed -g genome.txt > shufpromoter.bed
bedtools jaccard -a cpg.bed -b shufpromoter.bed



10. Which hESC ChromHMM state (e.g., 11_Weak_Txn, 10_Txn_Elongation) represents the most number of base pairs in the genome? (Hint: you will need to use awk or perl here, as well as the groupby tool.)


bedtools groupby -g $4 -c $

awk 'BEGIN { FS=OFS="\t" } {print $0,$3-$2}' hesc.chromHmm.bed | bedtools groupby -g 4 -c 5 -o sum 

