# Capture-C data processing
### CCanalyser3 pipeline: alignment and generating bedgraphs

Please refer to the [Capture-C paper: Hughes, J. R. et al. (2014)](https://www.nature.com/articles/ng.2871) and [its GitHub](https://github.com/Hughes-Genome-Group/captureC/releases) for details on the experimental technique and its alignment/analysis tool: CCanalyser.

Requires the folowing modules:
``` sh
PATH=$PATH:~/.local/bin
chmod u=rwx *.pl

perl -e module load trim_galore/0.3.1
perl -e module load bowtie/1.0.0
perl -e module load ucsctools/1.0
perl -e module load flash/1.2.8
```

Several additional files are required:
1. [CCanalyser3](https://upenn.box.com/s/yr7si2ulujrzos60t4v3msbw8wmt5snx)
2. [Chromosome sizes](https://upenn.box.com/s/2dbdvnzxwyx3m69o3lj1b281a5iiyajz)
3. [Genome fragments (DpnII)](https://upenn.box.com/s/ymdywbt0fygjcpcc33p51twrm3mmpcun)
4. [Bowtie indexed hg19](https://upenn.box.com/s/pdgvs4i6n1pnki7jnmpe41d5j1v4d2qk)
5. [Capture oligo coordinates](https://upenn.box.com/s/tmu1x9yzem9160et7bl1bqnhrx10yi54)

Use [WT, 5-loci Capture-C data as an example](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4077757)
``` sh
nohup trim_galore --paired 1887_Read1.fq 1887_Read2.fq; 
nohup ~/FLASH-1.2.8/flash --interleaved-output 1887_Read1_val_1.fq 1887_Read2_val_2.fq; 
cat out.notCombined.fastq out.extendedFrags.fastq > 1887_combined_reads.fastq; 

nohup ./dpnII2E.pl 1887_combined_reads.fastq; 

nohup bowtie -p 1 -m 2 --best --strata --sam --chunkmb 256 hg19 1887_combined_reads_REdig.fastq run1887_REdig.sam; 

perl CCanalyser3_mod.pl -b "<path to genome sizes>" -f "<path to run1887_REdig.sam>" -r "<path to DpnII genome fragments>" -genome hg19 -o "<path to oligo coordinates>" -s run_1887;
```
Once complete, bedgraphs can be extracted from .gff output files for downstream analysis. 


### Combine replicates and normalize

Take two WT replicates at C21S4 as an example. In R:
``` r
library('plyr')
wt_1 <- read.table("run1887_REdig_CC3_CapC_C21_S4.bedgraph", skip=1)
wt_2 <- read.table("run1888_REdig_CC3_CapC_C21_S4.bedgraph", skip=1)
wt_comb <- rbind (wt_1,wt_2)
wt_comb_norm <- ddply(wt, wt("V1","V2","V3"), numcolwise(sum))
wt_comb_norm$V4 <- wt_comb_norm$V4/sum(wt_comb_norm$V4)

write.table(wt_comb_norm, file="WT_combined_normed_C21_S4.bedgraph",quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
rm(list = ls())
``` 


### Compute Directionality Index across distance ranges

A folder should contain [Capture-C bedgraphs](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137373) of all genotypes for a genetically perturbed insertion locus. 

File needed to specify the insertion locus of interest: [C21S4, for example](https://upenn.box.com/s/a1hn241t1ncgqrco41g6zat36lzppt4c).

In R:
``` r
setwd("<path to working directory>")
file.names <- dir(getwd(), pattern=".bedgraph")

insert_loci_filename <- "insertion_loci_c21s4.txt"

DI_output <- matrix(, nrow = 0, ncol = 3)

for (i in 1:length(file.names)){
  bg <- read.table(file.names[i])
  f <- strsplit(file.names[i], "\\.")[[1]][1]
  
  insert_loci <- read.table(insert_loci_filename)
  insert_site <- insert_loci[1,][1]
  chrom <- insert_loci[1,][2]
  start_coord <- insert_loci[1,][3]
  end_coord <- insert_loci[1,][4]
  
  insert_chrom_frags <- bg[as.character(bg$V1) == as.character(chrom[1,1]),]
  for (j in c(50000, 100000, 250000, 500000, 1000000)) {
    upstream_int <- insert_chrom_frags[insert_chrom_frags$V2 >= (as.numeric(start_coord) - j) & insert_chrom_frags$V3 <= as.numeric(start_coord), ]
    downstream_int <- insert_chrom_frags[insert_chrom_frags$V2 >= as.numeric(end_coord) & insert_chrom_frags$V3 <= (as.numeric(end_coord) + j), ]
    
    downstream_sum <- sum(downstream_int$V4)
    upstream_sum <- sum(upstream_int$V4)
    
    avg_int <- (downstream_sum + upstream_sum)/2
    
    direct_idx <- (downstream_sum - upstream_sum)/abs(downstream_sum - upstream_sum) * ((upstream_sum - avg_int)^2/avg_int + (downstream_sum - avg_int)^2/avg_int)
    
    DI_output <- rbind(DI_output, c(f, j, direct_idx))
  }
}

f_name <- paste(as.character(insert_site[1,1]), "DI", "txt", sep=".")
write.table(DI_output, file=f_name, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
``` 
