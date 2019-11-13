# assumes you have SSH'ed and qlogin'ed
```
cd wed
cd mydata

slide 17
curl https://s3.amazonaws.com/gemini-tutorials/trio.trim.vep.vcf.gz > trio.trim.vep.vcf.gz
curl https://s3.amazonaws.com/gemini-tutorials/denovo.ped > denovo.ped
gemini load --cores 2 \
            -v  trio.trim.vep.vcf.gz \
            -t VEP \
           --tempdir . \
           --skip-gene-tables --skip-cadd --skip-gerp-bp \
            -p denovo.ped \
        trio.trim.vep.denovo.db

curl http://home.chpc.utah.edu/~u1138933/gemini_db/trio.trim.vep.denovo.db > trio.trim.vep.denovo.db
```

# slide 19
```
gemini de_novo trio.trim.vep.denovo.db
```

# type Ctrl+C to stop output if you'd like (should take 20 seconds to complete)

# slide 21
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" trio.trim.vep.denovo.db

# slide 23
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" trio.trim.vep.denovo.db | wc -l

# slide 25
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
trio.trim.vep.denovo.db | wc -l

# slide 26
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
--min-gq 20 \
trio.trim.vep.denovo.db | wc -l

# slide 28: Counting in not passed variants - means with warning due to paralogous genes, PCR eerrors
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
--min-gq 20 \
--filter "filter is NULL" \
trio.trim.vep.denovo.db | wc -l


gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" -d 6 --min-gq 20  --filter "filter is NULL and impact_severity != 'LOW' and max_aaf_all<0.0002941176471" trio.trim.vep.denovo.db
chrom	start	end	ref	alt	filter	qual	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count
chr17	76471351	76471352	G	A	None	1245.9699707	DNAH17	splice_region_variant	14284	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	G/G,G/G,G/A	4805	1
chr22	43027436	43027437	C	T	None	1320.0300293	CYB5R3	missense_variant	16718	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	C/C,C/C,C/T	4805	1
ubuntu@ip-172-31-4-4:~/wed/mydata$ gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" -d 6 --min-gq 20  --filter "filter is NULL and impact_severity != 'LOW' and max_aaf_all<0.005" trio.trim.vep.denovo.db
chrom	start	end	ref	alt	filter	qual	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count
chr17	76471351	76471352	G	A	None	1245.9699707	DNAH17	splice_region_variant	14284	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	G/G,G/G,G/A	4805	1
chr22	43027436	43027437	C	T	None	1320.0300293	CYB5R3	missense_variant	16718	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	C/C,C/C,C/T	4805	1
chr15	41229630	41229631	T	G	None	2116.48999023	DLL4	missense_variant	7892	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	T/T,T/T,T/G	4805	1
chr22	23040690	23040691	G	A	None	2574.36010742	IGLV2-23	missense_variant	15541	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	G/G,G/G,G/A	4805	1
chr22	23040728	23040729	G	C	None	2087.54003906	IGLV2-23	missense_variant	15542	family1	1805(1805;unaffected;female),1847(1847;unaffected;male),4805(4805;affected;male)	G/G,G/G,G/C	4805	1

# slide 30
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
--min-gq 20 \
--filter "(filter is NULL or filter=='SBFilter')" \
trio.trim.vep.denovo.db | wc -l

# slide 32
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
--min-gq 20 \
--filter "(filter is NULL or filter=='SBFilter') and impact_severity != 'LOW'" \
trio.trim.vep.denovo.db | wc -l

# slide 35
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
--min-gq 20 \
--filter "(filter is NULL or filter=='SBFilter') and impact_severity != 'LOW' and max_aaf_all <= 0.005" \
trio.trim.vep.denovo.db | wc -l

# slide 35
gemini de_novo --columns "chrom, start, end, ref, alt, filter, qual, gene, impact" \
-d 6 \
--min-gq 20 \
--filter "(filter is NULL or filter=='SBFilter') and impact_severity != 'LOW' and max_aaf_all <= 0.005" \
trio.trim.vep.denovo.db
