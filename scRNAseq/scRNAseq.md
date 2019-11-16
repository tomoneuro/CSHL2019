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
# Post-cellranger analysis and interpretation

# Seurat V3 workflow overview
for a single sample or pre-integrated samples

```
scrna.counts <- Read10X(data.dir = "/yourpath/outs/filtered_feature_bc_matrix") scrna <- CreateSeuratObject(counts = scrna.counts)
```

IF, you eliminate genes with essentially no expression, and cells with very few genes,
Caution!
• Rare genes may be expressed in only a few cells
• Different cell types have different numbers of genes. Example: Red blood cells express
only ~200 genes in some 10x data sets. -> Don’t filter out potentially important cells.
```
scrna <- CreateSeuratObject(counts = scrna.counts, min.cells = 10, min.features = 100, project = Project1)
```
```
scrna <- SCTransform(scrna) # Find variable genes, scale, and normalize
scrna <- NormalizeData(object = scrna) # older pipeline
scrna <- FindVariableFeatures(object = scrna) # older pipeline scrna <- ScaleData(object = scrna) # older pipeline
scrna <- RunPCA(object = scrna) # Principal Component Analysis scrna <- FindNeighbors(object = scrna)
scrna <- FindClusters(object = scrna)
scrna <- RunTSNE(object = scrna)
scrna <- RunUMAP(object = scrna) DimPlot(object = scrna, reduction = "tsne") UMAPPlot(object = scrna)
```

# Calculate percentage of mitochondrial genes
```
mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna), value = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'));
scrna[['percent.mito']] <- percent.mito; # assign it to the meta data
```
# ribosomal genes
```
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna), value = TRUE); percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'));
scrna[['percent.ribo']] <- percent.ribo; # assign it to the meta data
vln <- VlnPlot(object = scrna, features = c("percent.mito", "percent.ribo"), pt.size=0, ncol = 2, group.by="Batch"); # make a violin plot, and color by batch or sample
vln <- VlnPlot(object = scrna, features = "nCount_RNA", pt.size=0, group.by="Batch", y.max=25000) vln <- VlnPlot(object = scrna, features = "nFeature_RNA", pt.size=0, group.by="Batch")
```

# Filter data
Remove cells with high mitochondrial content, too many genes, too few genes, etc
```
scrna <- subset(x = scrna, subset = nFeature_RNA > 200 & nFeature_RNA < Feature95 & percent.mito < mc.hi)
```

# Calculate cell cycle phase of each cell
Cell cycle signature can dominate tSNE/UMAP plots and may need to be removed.
Use the plots on the previous slides to choose appropriate thresholds for your data
• Some example thresholds:
• nFeature_RNA between 200 and 4000, or between 200 and 95th percentile (for example)
• nFeature > x, nUMI < 93rd percentile
• Mitochondrial content below 5% or 10%
• Ribosomal content below 50%

```
cell.cycle.tirosh <- read.table("CellCycleTirosh.txt", sep='\t', header=FALSE); s.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G1/S")]; g2m.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G2/M")];
scrna <- CellCycleScoring(object=scrna, s.features=s.genes, g2m.features=g2m.gen es, set.ident=FALSE)
```

# Calculate cell cycle phase of each cell
Cell cycle signature can dominate tSNE/UMAP plots and may need to be removed. Seurat has a built-in function that calculates relative expression level of G1/S and G2/M genes defined in Tirosh et al 2016.
```
cell.cycle.tirosh <- read.table("CellCycleTirosh.txt", sep='\t', header=FALSE); s.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G1/S")]; g2m.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G2/M")];
scrna <- CellCycleScoring(object=scrna, s.features=s.genes, g2m.features=g2m.gen es, set.ident=FALSE)
```

# Step 5: Normalize, scale, control for unwanted variation
Goal: remove technical effects while preserving biological variation
Older normalization approach(es) scale every gene in the cell by the same factor
(sequencing depth)
• NormalizeData:
	• Divide by total counts in each cell
	• Scale to fixed counts (default is 1x104 use 1x106 for CPM)
	• Add1
	• natural log
• FindVariableFeatures: use a variance stabilizing transformation
• ScaleData: subtract mean, divide by standard deviation, remove unwanted signal using multiple regression

## SCTransform Normalization/Scaling function
```
# SCTransform (with removal of cell cycle signal):
scrna <- SCTransform(scrna, vars.to.regress = c("S.Score", "G2M.Score"), verbose=FALSE) # Note that this replaces the three separate steps implemented in older workflows:
scrna <- NormalizeData(object = scrna, normalization.method = "LogNormalize", scale.factor = 1e4) # fea ture counts divided by total, multiplied by scale factor, add 1, ln-transformed
scrna <- FindVariableFeatures(object = scrna, selection.method = 'vst', mean.cutoff = c(0.1,8), dispers ion.cutoff = c(1, Inf)) # designed to find ~2000 variable genes
scrna <- ScaleData(object = scrna, features = rownames(x = scrna), vars.to.regress = c("S.Score","G2M.S core"), display.progress=FALSE) # center and regress out unwanted variation
```

For class project, SeqTech2019_Data, start with following code to open the file.
```
scrna.counts <- Read10X_h5("./Desktop/SeqTech2019_Data/filtered_feature_bc_matrix.h5")

cat test.fa
>scaff1
AAAAAAAAAA
>scaff2
TTTTTTTTT
1:08 PM
samtools faidx test.fa
1:08 PM
cat test.fa.fai
scaff1    10    8    10    11
scaff2    9    27    9    10
```
