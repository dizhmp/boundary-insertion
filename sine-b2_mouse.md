# SINE B2: possible context dependency in forming putative new domains

## Prepare peak sets, boundaries, and TSSs. 

###[SINE B2-associated CTCF peaks](https://upenn.box.com/s/bgd8c0dsebg92u8lcsw81e9pnd81ouun). This was generated through the following steps:

1. [Thybert et al, 2018](https://genome.cshlp.org/content/28/4/448.long)
2. Download [Mus musculus repeat-associated CTCF peak set](https://www.ebi.ac.uk/research/flicek/publications/FOG21)
3. Filter down to B2 RepeatFamily, and convert to BED format. 
4. Use [UCSC LiftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert from mm10 to mm9.


###[5-way conserved (ancestral) CTCF peaks](https://upenn.box.com/s/q1ze0a1ka7fht4koohnn0wa7ezbyyyjq). This was derived through the following steps:

1. [Schmidt et al, 2012](https://www.cell.com/fulltext/S0092-8674(11)01507-8)
2. Download [5-way (placental) shared CTCF binding events (mm9)] http://ftp.ebi.ac.uk/pub/databases/vertebrategenomics/FOG03/
3. Convert it to BED format.


###[All CTCF peaks (erythroid)](https://upenn.box.com/s/iyocw3vyt9dsscmo0iu2iuaavdbyd0z7). This was obtained through the following steps: 

1. [Zhang et al, 2019](https://www.nature.com/articles/s41586-019-1778-y)
2. Download Supplementary Table 2 under Supplementary Information
3. Filter by "CTCF_reps_merge_Asyn = 1" to obtain all CTCF peaks under asynchronized cell cycle condition.
4. Convert it to BED format. 


###[Domain boundaries (erythroid)](https://upenn.box.com/s/e9b32le0jvsv8imuv2fnprlqs0fg1dv3), obtained through the following steps:

1. [Zhang et al, 2019](https://www.nature.com/articles/s41586-019-1778-y)
2. Download Supplementary Table 1 under Supplementary Information
3. Use "Final_emergent_boundaries," and extend 10 kb upstream and downstream.
4. Convert it to BED format.


###[TSS annotation](https://upenn.box.com/s/rgnnzgrrs40epicvtbmhh3xuo7ras4ft), obtained through the following steps:

1. [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables)
2. Select the following set of parameters:

``` sh
	# http://genome.ucsc.edu/cgi-bin/hgTables
	# genome: mouse
	# assembly: mm9
	# group: gene and gene predictions
	# track: RefSeq genes
	# table: refgene
	# region: genome
	# intersection with knownGene: track: UCSC genes (knownGene); 80% overlap
	# get output: upstream by 1 base
	# saved as: refseq_ucsc_mm9.bed
```

3. Remove redundant coordinates and clean up the format:
``` sh
awk '!seen[$1,$2,$3]++' refseq_ucsc_mm9.bed > refseq_ucsc_mm9_unique.bed
bedtools sort -i refseq_ucsc_mm9_unique.bed > refseq_ucsc_mm9_unique_sorted.bed
```

Further refine the TSS annotation:
In R:

``` r
library('dplyr')
refseq_ucsc_mm9_unique_sorted <- read.table("refseq_ucsc_mm9_unique_sorted.bed")

# extract refseq ids
for (i in 1:nrow(refseq_ucsc_mm9_unique_sorted)) {
  refseq_ucsc_mm9_unique_sorted[i,]$V5 <- strsplit(toString(refseq_ucsc_mm9_unique_sorted[i,]$V4), "_up_")[[1]][1] 
}

refseq_ucsc_mm9_unique_sorted$V4 <- NULL

# unique refseq:
refseq_count <- as.data.frame(table(refseq_ucsc_mm9_unique_sorted$V5))
colnames(refseq_count)[1] <- "V5"

# u: unique; s: strand; occu: occurance.
refseq_uc_mm9_u_s_occu <- merge(refseq_ucsc_mm9_unique_sorted, refseq_count, by="V5", sort = FALSE)[, union(names(refseq_ucsc_mm9_unique_sorted), names(refseq_count))]

refseq_uc_mm9_u_s_redun <- refseq_uc_mm9_u_s_occu[refseq_uc_mm9_u_s_occu$Freq>1, ]

redun_refseq <- unique(refseq_uc_mm9_u_s_redun$V5)
dis_refseq_redun <- c()
redunt_refseq_set <- refseq_uc_mm9_u_s_redun[refseq_uc_mm9_u_s_redun$V5==toString(redun_refseq[1]),]
uniq_from_redun_refseq <- redunt_refseq_set[0,]

for (i in 1:length(redun_refseq)) {
  redunt_refseq_set <- refseq_uc_mm9_u_s_redun[refseq_uc_mm9_u_s_redun$V5==toString(redun_refseq[i]),]
  dis_refseq_redun <- append(dis_refseq_redun, abs(min(redunt_refseq_set$V2) - max(redunt_refseq_set$V2)))
  # check if redudant refseq ids map to the same chromosome (YES --> check distance; No --> include all)
  if (all.equal(redunt_refseq_set$V1, rev(redunt_refseq_set$V1)) != TRUE){
    uniq_from_redun_refseq <- rbind(uniq_from_redun_refseq, redunt_refseq_set)}
  else {
    # if on different strands, keep all
    if (all.equal(redunt_refseq_set$V6, rev(redunt_refseq_set$V6)) != TRUE) {
      uniq_from_redun_refseq <- rbind(uniq_from_redun_refseq, redunt_refseq_set)
    }
    else {
    # if same strands, save the upstream one
    if (redunt_refseq_set$V6[1] == "+") {
      uniq_from_redun_refseq <- rbind(uniq_from_redun_refseq, redunt_refseq_set[redunt_refseq_set$V2 == min(redunt_refseq_set$V2),])}
    else {
      uniq_from_redun_refseq <- rbind(uniq_from_redun_refseq, redunt_refseq_set[redunt_refseq_set$V2 == max(redunt_refseq_set$V2),])
    }
    }
  }
}

refseq_ucsc_mm9_refseq <- rbind(refseq_uc_mm9_u_s_occu[refseq_uc_mm9_u_s_occu$Freq == 1,], uniq_from_redun_refseq)

refseq_ucsc_mm9_refseq$dist_from_next <- 0
refseq_ucsc_mm9_refseq$remove <- "N"

# only keep 1 tss on the same strand within a 100bp window.
for (i in 2:nrow(refseq_ucsc_mm9_refseq)) {
  dis <- refseq_ucsc_mm9_refseq[i,]$V2 - refseq_ucsc_mm9_refseq[i-1,]$V2
  if ((dis > 0) && (dis <= 100)) {
    if (refseq_ucsc_mm9_refseq[i,]$V6 == refseq_ucsc_mm9_refseq[i-1,]$V6) {
      refseq_ucsc_mm9_refseq[i,]$remove <- "Y"
    }
  }
  refseq_ucsc_mm9_refseq[i-1,]$dist_from_next <- dis
}

refseq_ucsc_mm9_refseq <- refseq_ucsc_mm9_refseq[refseq_ucsc_mm9_refseq$remove == "N",]

colnames(refseq_ucsc_mm9_refseq)[1] <- "chr"
colnames(refseq_ucsc_mm9_refseq)[2] <- "start"
colnames(refseq_ucsc_mm9_refseq)[3] <- "end"
colnames(refseq_ucsc_mm9_refseq)[4] <- "name"

refseq_ucsc_mm9_refseq_bed4 <- refseq_ucsc_mm9_refseq[, 1:4]
write.table(refseq_ucsc_mm9_refseq_bed4, file = "refseq_ucsc_refined_relaxed.bed", quote = FALSE, sep= "\t", row.names = FALSE, col.names = TRUE)
```

remove carrige returns:
``` sh
sed -i "s/^M//g" refseq_ucsc_refined_relaxed.bed
```

## Overlap analysis

SINE B2-associated CTCF peaks present in at least two tissues:
``` sh
bedtools intersect -a mus_B2_mm9_odom_flicek_unique.bed -b JC4_CTCF_InAsyn.bed -u -wa > CTCF_B2_liver_JC4_FinalAsyn_overlap.bed
```

Ancestral CTCF peaks present in at least two tissues:
``` sh
bedtools intersect -a CTCF_mm9_5way.bed -b JC4_CTCF_InAsyn.bed -u -wa > CTCF_5way_liver_JC4_FinalAsyn_overlap.bed
```

Boundaries with ancestral CTCF peaks:
``` sh
bedtools intersect -a boundaries_final.bed -b CTCF_5way_liver_JC4_FinalAsyn_overlap.bed -u -wa > JC4_boundaries_5wayliverJC4_FinalAsyn.bed
```

Boundaries with SINE B2 CTCF peaks and without ancestral CTCF:
``` sh
bedtools intersect -a boundaries_final.bed -b CTCF_B2_liver_JC4_FinalAsyn_overlap.bed -u -wa > JC4_boundaries_B2liverJC4_FinalAsyn.bed

bedtools intersect -a JC4_boundaries_B2liverJC4_FinalAsyn.bed -b JC4_boundaries_5wayliverJC4_FinalAsyn.bed -v -wa > JC4_boundaries_B2liverJC4_Exclude5way_FinalAsyn.bed
```

Extend a 200 kb window around each CTCF peak, in R:
``` r
B2_liver_JC4 <- read.table('CTCF_B2_liver_JC4_FinalAsyn_overlap.bed')
B2_liver_JC4$V4 <- rowMeans(B2_liver_JC4[,2:3])

B2_liver_JC4$V5 <- floor(B2_liver_JC4$V4 - 100000)
B2_liver_JC4$V6 <- ceiling(B2_liver_JC4$V4 + 100000)

B2_liver_JC4_200k <- B2_liver_JC4[,c(1,5,6)]
write.table(B2_liver_JC4_200k, file='B2_liver_JC4_200k_03162020.bed', row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
```

Remove carrige returns:
``` sh
sed -i "s/^M//g" B2_liver_JC4_200k_03162020.bed
```

Count the # of TSSs within a 200 kb window for each SINE B2-associated CTCF peak:
``` sh
bedtools intersect -a B2_liver_JC4_200k_03162020.bed -b refseq_ucsc_refined_relaxed.bed -c -wa > B2_liver_JC4_200k_rel_TSS_03162020.bed

paste CTCF_B2_liver_JC4_FinalAsyn_overlap.bed B2_liver_JC4_200k_rel_TSS_03162020.bed | cut -f 1,2,3,8 > CTCF_B2_liver_JC4_TSSCounts_FinalAsyn.bed
```

Co-localization between SINE B2-associated CTCF peaks and boundaries with SINE B2 CTCF peaks and without ancestral CTCF:
``` sh
# need to change file names
bedtools intersect -a JC4_boundaries_B2liverJC4_Exclude5way_FinalAsyn.bed -b CTCF_B2_liver_JC4_TSSCounts_FinalAsyn.bed -c > bounds_counts_b2_ctcf.bed

bedtools intersect -a bounds_counts_b2_ctcf.bed -b CTCF_B2_liver_JC4_TSSCounts_FinalAsyn.bed -loj > bounds_counts_b2_ctcf_loj.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $5, $6, $7, $8, $4}' bounds_counts_b2_ctcf_loj.bed > bounds_counts_b2_ctcf_loj_rearr.bed
```

Co-localization by TSS density, in R:
``` r
library('plyr')

CTCF_B2_TSS <- read.table("CTCF_B2_liver_JC4_TSSCounts_FinalAsyn.bed")

colnames(CTCF_B2_TSS)[1] <- "chr"
colnames(CTCF_B2_TSS)[2] <- "start"
colnames(CTCF_B2_TSS)[3] <- "end"
colnames(CTCF_B2_TSS)[4] <- "rel_200k"
CTCF_B2_TSS$redun <- 0
CTCF_B2_TSS$in_boundary <- 0

CTCF_B2_Bounds <- read.table("bounds_counts_b2_ctcf_loj_rearr.bed")
colnames(CTCF_B2_Bounds)[1] <- "chr"
colnames(CTCF_B2_Bounds)[2] <- "start"
colnames(CTCF_B2_Bounds)[3] <- "end"
colnames(CTCF_B2_Bounds)[4] <- "rel_200k"
colnames(CTCF_B2_Bounds)[5] <- "redun"
CTCF_B2_Bounds$redun <- CTCF_B2_Bounds$redun - 1
CTCF_B2_Bounds[, 4] <- 0
CTCF_B2_Bounds$in_boundary <- 1

# count multiple b2 ctcfs in a single boundary
for (k in 1:nrow(CTCF_B2_Bounds)) {
  next_redun_rows <- CTCF_B2_Bounds[k,]$redun
  if (next_redun_rows > 0) {
    for (m in 1:next_redun_rows) {
      CTCF_B2_Bounds[k+m,]$redun <- 0
    }
  }
}

ctcf_b2_tss_bounds_pre <- rbind(CTCF_B2_TSS, CTCF_B2_Bounds)
ctcf_b2_tss_bounds <- ddply(ctcf_b2_tss_bounds_pre, colnames(ctcf_b2_tss_bounds_pre)[1:3], numcolwise(sum))

for (i in 4:(ncol(ctcf_b2_tss_bounds)-2)) {
  dist_tss <- colnames(ctcf_b2_tss_bounds)[i]
  tss_density_subset <- as.data.frame(table(ctcf_b2_tss_bounds[,i]))
  tss_density_subset$in_boundary <- 0
  tss_density_subset$redun <- 0
  for (j in 1:nrow(tss_density_subset)) {
    tss_den <- tss_density_subset[j,]$Var1
    ss <- ctcf_b2_tss_bounds[ctcf_b2_tss_bounds[,i]==tss_den,]
    tss_density_subset$in_boundary[j] <- sum(ss$in_boundary)
    tss_density_subset$redun[j] <- sum(ss$redun)
  }
  write.table(tss_density_subset, file = paste(dist_tss, "_IncluRedun_Asyn_0316.bed", sep=""), quote = FALSE, sep= "\t", row.names = FALSE, col.names = TRUE)
}
```
