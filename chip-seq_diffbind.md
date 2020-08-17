# ChIP-seq: CTCF/RAD21 differential binding 
### [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html):

After installing DiffBind in R following the link above, several files are needed for comparing Clone 21 vs. non-Clone 21 CTCF/RAD21 binding as an example:

1. Sample sheets [CTCF](https://upenn.box.com/s/89l13t9w7of44f1hjf24w3bkxulrckoy) and [RAD21](https://upenn.box.com/s/r67mmdfv0exwxg6vaocitnw3n512u6dn): paths to "bamReads" and "Peaks" need to be updated once these files have been downloaded. 
2. [Bam and peak files](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152721)

``` r
# need to modify the path
setwd("<path to sample sheets>")

library("DiffBind")

ctcf <- dba(sampleSheet = "ctcf_sample_sheet.csv")
rad21 <- dba(sampleSheet = "rad21_sample_sheet.csv")

write.table( x = dba.show(ctcf), file = "ctcf_import.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = dba.show(rad21), file = "rad21_import.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

# peak centering
ctcf <- dba.count(ctcf, minOverlap=2, summits=250)
rad21 <- dba.count(rad21, minOverlap=2, summits=250)

write.table( x = dba.show(ctcf), file = "ctcf_frip.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = dba.show(rad21), file = "rad21_frip.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

# a peak eligible for differential binding analysis should be identified in 2 samples
ctcf_1 <- dba.contrast(ctcf, minMembers=2, categories = DBA_CONDITION)
rad21_1 <- dba.contrast(rad21, minMembers=2, categories = DBA_CONDITION)

ctcf_2 <- dba.analyze(ctcf_1, bFullLibrarySize =  FALSE, bTagwise = FALSE)
rad21_2 <- dba.analyze(rad21_1, bFullLibrarySize =  FALSE, bTagwise = FALSE)

ctcf.DB <- dba.report(ctcf_2, contrast = c(1:1), bDB = TRUE)
rad21.DB <- dba.report(rad21_2, contrast = c(1:1), bDB = TRUE)

write.table( x = dba.show(ctcf.DB), file = "ctcf_report.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = dba.show(rad21.DB), file = "rad21_report.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

# becomes more useful for multiple pairwise comparisons 
for (i in 1: length(ctcf_2$contrasts)){
  ctcf_peaks.DB <- dba.report(ctcf_2, contrast = i)
  rad21_peaks.DB <- dba.report(rad21_2, contrast = i)
  
  ctcfDB_df = as(ctcf_peaks.DB, "data.frame")
  rad21DB_df = as(rad21_peaks.DB, "data.frame")
  
  write.table( x = ctcfDB_df, file = paste ("ctcf_db_", i, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
  write.table( x = rad21DB_df, file = paste ("rad21_db_", i, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

  # report all peaks and their quantifications
  ctcf_peaks.allpeaks <- dba.report(ctcf_2, contrast = i, th=1, bCounts=TRUE)
  rad21_peaks.allpeaks <- dba.report(rad21_2, contrast = i, th=1, bCounts=TRUE)
  
  ctcfALLPEAKS_df = as(ctcf_peaks.allpeaks, "data.frame")
  rad21ALLPEAKS_df = as(rad21_peaks.allpeaks, "data.frame")
  
  write.table( x = ctcfALLPEAKS_df, file = paste ("ctcf_allpeaks_", i, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
  write.table( x = rad21ALLPEAKS_df, file = paste ("rad21_allpeaks_", i, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
}
```


### Extract insertion-proximal peaks
Convert .csv to .bed in R:

``` r
file.names <- dir(getwd(), pattern=".csv")

for (i in 1:length(file.names)){
  csv <- read.table(file.names[i], sep = ",")
  f <- strsplit(file.names[i], "\\.")[[1]][1]
  f_name <- paste(f, "bed", sep=".")
  write.table(csv[-1,c(1:3, 9, 11)], file=f_name, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}
```


Use [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).
Requires distance ranges:
1. [Clone21 insertion-proximal range](https://upenn.box.com/s/owpgpaicrh76pa65r4q8aivbgkrwjrsf)
2. [Clone25 insertion-proximal range](https://upenn.box.com/s/06eyyzwctou2m9h1yl3yj42zo83t8b36)
``` sh
for f in "$PWD"/*.bed
do
	bedtools intersect -wa -a "$f" -b C21_3Mb.txt > "$f.C21_3Mb.bed"
	bedtools intersect -wa -a "$f" -b C25_3Mb.txt > "$f.C25_3Mb.bed"
done
```


Obtain full quantifications for insertional-proximal peaks in R:
``` r
file.names <- dir(getwd(), pattern="allpeaks_1.csv")

for (i in 1:length(file.names)){
  csv <- read.table(file.names[i], sep = ",", header = TRUE)
  f <- strsplit(file.names[i], "\\.")[[1]][1]
  f_name <- paste(f, "bed.C21_3Mb.bed", sep=".")
  c21_3m <- read.table(f_name, sep='\t')
  colnames(c21_3m)[1] <- 'seqnames'
  colnames(c21_3m)[2] <- 'start'
  colnames(c21_3m)[3] <- 'end'
  colnames(c21_3m)[4] <- 'Fold'
  colnames(c21_3m)[5] <- 'FDR'
  csv_c21_3m <- merge(csv, c21_3m)
  f_out <- paste(f, "C21_3Mb.tsv", sep="_")
  write.table(csv_c21_3m, file=f_out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}
```
