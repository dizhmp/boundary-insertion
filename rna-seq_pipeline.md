# RNA-seq: differential expression analysis

### Transcript-level quantification using [Salmon](https://www.nature.com/articles/nmeth.4197)

[Salmon documentation](https://combine-lab.github.io/salmon/getting_started/) and [releases](https://github.com/COMBINE-lab/salmon/releases). We used v0.9.1. for this analysis. 

After installation, the following items are required:
1. Reference transcriptome ([Ensembl hg19: Homo_sapiens.GRCh37.67.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-67/fasta/homo_sapiens/cdna/))
2. [RNA-seq raw reads](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137374)

Index the reference transcriptome:
``` sh
# indexing hg19 transcriptome
./Salmon-latest_linux_x86_64/bin/salmon index -t hg19_transcriptome.tar.gz -i hg19_transcriptome_index
```

Transcript-level quantification, using WT, Clone21 and Clone25 as an example:
``` sh
# salmon quant rna seq
./Salmon-latest_linux_x86_64/bin/salmon quant -i hg19_transcriptome_index -l A -1 1732_Read1.fq.gz -2 1732_Read2.fq.gz -p 8 -o quants/1732;
./Salmon-latest_linux_x86_64/bin/salmon quant -i hg19_transcriptome_index -l A -1 1733_Read1.fq.gz -2 1733_Read2.fq.gz -p 8 -o quants/1733; 
./Salmon-latest_linux_x86_64/bin/salmon quant -i hg19_transcriptome_index -l A -1 1734_Read1.fq.gz -2 1734_Read2.fq.gz -p 8 -o quants/1734; 
./Salmon-latest_linux_x86_64/bin/salmon quant -i hg19_transcriptome_index -l A -1 1735_Read1.fq.gz -2 1735_Read2.fq.gz -p 8 -o quants/1735; 
./Salmon-latest_linux_x86_64/bin/salmon quant -i hg19_transcriptome_index -l A -1 1736_Read1.fq.gz -2 1736_Read2.fq.gz -p 8 -o quants/1736; 
./Salmon-latest_linux_x86_64/bin/salmon quant -i hg19_transcriptome_index -l A -1 1737_Read1.fq.gz -2 1737_Read2.fq.gz -p 8 -o quants/1737;
```


### Differential expression analysis using [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

Requires:
1. [Sample sheet: samples.txt](https://upenn.box.com/s/vkm5a6754ukgqdy0y047ggyqast7hcbj)


In R, import transcript-level estimates towards gene-level quantification using [tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html): 
``` r
# modify working directory
setwd("<path to "quants" generated using Salmon above>")

# install and load required modules
library("tximport")
library("readr")

# read sample sheet
samples <- read.table("samples.txt", header=TRUE)

# compare Clone 21 vs. non-Clone 21 
samples$condition <- factor(rep(c("NoC21Insertions","C21Insertions","NoC21Insertions"),each=2))
rownames(samples) <- samples$run
files <- file.path(samples$run, "quant.sf")
names(files) <- samples$run

# specify transcriptome
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
txdf <- transcripts(edb, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
```

Run DESeq2:
``` r
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

# pre-filtering: only keep rows with at least 10 reads
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]

# specify reference condition:
dds$condition <- relevel(dds$condition, ref = "NoC21Insertions")

# run DESeq
dds <- DESeq(dds)

# running default parameters
res <- results(dds)

# results summary
summary(res)
```


### Obtain quantifications on insertion-proximal genes

Requires:
1. [Distance ranges from insertions](https://upenn.box.com/s/dnbm7bm2dted0tt0ebprox4coqm4x0rj)


```r
# adding gene names
library(Homo.sapiens)
library("biomaRt")
res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )


ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl )

idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
res$entrez <- genemap$entrezgene[idx]

# made .txt files with certain window sizes near insertions
# make GRanges subjucts from .txt
C21_1bp <- read.table("C21_1bp.txt")
as.data.frame(C21_1bp)
colnames(C21_1bp)[1] <- "chr"
colnames(C21_1bp)[2] <- "start"
colnames(C21_1bp)[3] <- "end"
C21_1bp.gr <- makeGRangesFromDataFrame(C21_1bp)

C21_100kb <- read.table("C21_100kb.txt")
as.data.frame(C21_100kb)
colnames(C21_100kb)[1] <- "chr"
colnames(C21_100kb)[2] <- "start"
colnames(C21_100kb)[3] <- "end"
C21_100kb.gr <- makeGRangesFromDataFrame(C21_100kb)

C21_1Mb <- read.table("C21_1Mb.txt")
as.data.frame(C21_1Mb)
colnames(C21_1Mb)[1] <- "chr"
colnames(C21_1Mb)[2] <- "start"
colnames(C21_1Mb)[3] <- "end"
C21_1Mb.gr <- makeGRangesFromDataFrame(C21_1Mb)

C21_3Mb <- read.table("C21_3Mb.txt")
as.data.frame(C21_3Mb)
colnames(C21_3Mb)[1] <- "chr"
colnames(C21_3Mb)[2] <- "start"
colnames(C21_3Mb)[3] <- "end"
C21_3Mb.gr <- makeGRangesFromDataFrame(C21_3Mb)

C21_10Mb <- read.table("C21_10Mb.txt")
as.data.frame(C21_10Mb)
colnames(C21_10Mb)[1] <- "chr"
colnames(C21_10Mb)[2] <- "start"
colnames(C21_10Mb)[3] <- "end"
C21_10Mb.gr <- makeGRangesFromDataFrame(C21_10Mb)

# overlap GRanges with genes
C21_1bp_subset <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), C21_1bp.gr)
C21_100kb_subset <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), C21_100kb.gr)
C21_1Mb_subset <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), C21_1Mb.gr)
C21_3Mb_subset <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), C21_3Mb.gr)
C21_10Mb_subset <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), C21_10Mb.gr)

C21_100kb_subset_unique <- C21_100kb_subset[-match(C21_1bp_subset, C21_100kb_subset),]
C21_1Mb_subset_unique <- C21_1Mb_subset[-match(C21_100kb_subset, C21_1Mb_subset),]
C21_3Mb_subset_unique <- C21_3Mb_subset[-match(C21_1Mb_subset, C21_3Mb_subset),]
C21_10Mb_subset_unique <- C21_10Mb_subset[-match(C21_3Mb_subset, C21_10Mb_subset),]

C21_1bp_subset_in_DESeqResults <- match(res$entrez, C21_1bp_subset $ gene_id)
C21_100kb_subset_in_DESeqResults <- match(res$entrez, C21_100kb_subset_unique $ gene_id)
C21_1Mb_subset_in_DESeqResults <- match(res$entrez, C21_1Mb_subset_unique $ gene_id)
C21_3Mb_subset_in_DESeqResults <- match(res$entrez, C21_3Mb_subset_unique $ gene_id)
C21_10Mb_subset_in_DESeqResults <- match(res$entrez, C21_10Mb_subset_unique $ gene_id)

res $ within_C21_InsWin1bp <- C21_1bp_subset $ gene_id[C21_1bp_subset_in_DESeqResults] !="NA"
res $ within_C21_InsWin100kb <- C21_100kb_subset_unique $ gene_id[C21_100kb_subset_in_DESeqResults] !="NA"
res $ within_C21_InsWin1Mb <- C21_1Mb_subset_unique $ gene_id[C21_1Mb_subset_in_DESeqResults]  !="NA"
res $ within_C21_InsWin3Mb <- C21_3Mb_subset_unique $ gene_id[C21_3Mb_subset_in_DESeqResults]  !="NA"
res $ within_C21_InsWin10Mb <- C21_10Mb_subset_unique $ gene_id[C21_10Mb_subset_in_DESeqResults] !="NA"

# store row names/ensembl gene ids
genes_ensemblID_1bp_within_inserts <- rownames(res) [!is.na(res $ within_C21_InsWin1bp)]
genes_ensemblID_100kb_within_inserts <- rownames(res) [!is.na(res $ within_C21_InsWin100kb)]
genes_ensemblID_1Mb_within_inserts <- rownames(res) [!is.na(res $ within_C21_InsWin1Mb)]
genes_ensemblID_3Mb_within_inserts <- rownames(res) [!is.na(res $ within_C21_InsWin3Mb)]
genes_ensemblID_10Mb_within_inserts <- rownames(res) [!is.na(res $ within_C21_InsWin10Mb)]

# MA-plot
# circle genes
plotMA(res, ylim=c(-10,10), cex=0.75, alpha=0.009999, main = "Transcriptome Changes in Insertional Clone 21 (Red Dots Denoting Differentially Expressed Genes at FDR < 0.01)")
abline(h=c(-1,1), col="orange", lwd=2)

with(res[genes_ensemblID_1bp_within_inserts, ], {
  points(baseMean, log2FoldChange, col=rainbow(5)[5], cex=2, lwd=2)
})

with(res[genes_ensemblID_100kb_within_inserts, ], {
  points(baseMean, log2FoldChange, col=rainbow(5)[2], cex=2, lwd=2)
})

with(res[genes_ensemblID_1Mb_within_inserts, ], {
  points(baseMean, log2FoldChange, col=rainbow(5)[3], cex=2, lwd=2)
})

with(res[genes_ensemblID_3Mb_within_inserts, ], 
  {points(baseMean, log2FoldChange, col=rainbow(5)[4], cex=2, lwd=2)
})

legend("topright", c("Gene Body","<50Kb", "50kb-500kb", "500kb-1.5Mb"), fill=rainbow(5)[c(5,2,3,4)], horiz=FALSE, cex=0.7, y.intersp=0.6, text.width=0.8)

# write results to csv
write.csv(as.data.frame(res), file="C21_DESeq_res_apr2020_rerun_grch37.csv")
```
