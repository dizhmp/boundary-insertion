# Hi-C data processing
### HiC-Pro pipeline

Perhaps the easiest way to run [HiC-Pro](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x) is to use a [Singularity](https://singularity.lbl.gov/) image.

After installing [Singularity](https://singularity.lbl.gov/), download a [HiC-pro image (2.11.3-beta)](https://upenn.box.com/s/nfdhd5ne0qyxiagsnfex8cshb54c059f), or [a more recent version](https://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/).

Next, download (and extract if necessary) the following required files:

1. [Bowtie 2 indexed hg19](https://upenn.box.com/s/yv66sm6son68ot31wuvb12ctimdlzh9r)
2. [Genome sizes](https://upenn.box.com/s/elm36k2js9i5k5503oh0e78t5zu5dlxk)
3. [Genome fargments (DpnII)](https://upenn.box.com/s/bsp6gy3ojtprnp4vth8q6wwfmbtr3wxd)
4. [Config file](https://upenn.box.com/s/r3r74rdle1w0xi3jjvo7s8ggebmwzwi6)
5. [Corresponding Hi-C reads](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137376)

Open the config file, and modify lines 16, 20, 40, 52, 66 accordingly.

Please place raw sequencing reads for each sample/genotype into a separate subfolder, as required by HiC-Pro.

	singularity exec --bind <working directory path> hicpro_devel_ubuntu.img HiC-Pro -i <path to raw reads> -o <output directory> -c <path to the config file (e.g. config_03262020_crispr_full_git.txt)>

### Global contact matrix comparison

In R, Using [WT and Clone 25 matrices](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137372) as an example:

``` r
library(dplyr)
a <- read.table('hap1_wt_hic_20000_iced.matrix')
b <- read.table('hap1_clone25_hic_20000_iced.matrix')
c <- inner_join(a,b,by=c("V1","V2"))
```

Normalize/scaling by division by (# of valid interactions/approx. genomic bins)

Compute the Pearson correlation coefficient

``` r
c$V3.x <- c$V3.x/(169604428/154770)
c$V3.y <- c$V3.y/(169069685/154770)
cor(c$V3.x,c$V3.y)
```

Generate smooth scatter plots using log10 normalized contacts
``` r
c$V3.x <- log10(c$V3.x)
c$V3.y <- log10(c$V3.y)
smoothScatter(c$V3.x, c$V3.y, xlim = c(-5,0), ylim = c(-5,0))
```

### Heatmap generation

Extract Clone 21 insertion-proximal contact matrices, for example:

``` sh
awk '($1 >= 7857) && ($1 <= 8157) && ($2 >= 7857) && ($2 <= 8157)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site1.matrix
awk '($1 >= 7857) && ($1 <= 8157) && ($2 >= 7857) && ($2 <= 8157)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site1.matrix
awk '($1 >= 7857) && ($1 <= 8157) && ($2 >= 7857) && ($2 <= 8157)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site1.matrix

awk '($1 >= 12986) && ($1 <= 13286) && ($2 >= 12986) && ($2 <= 13286)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site2.matrix
awk '($1 >= 12986) && ($1 <= 13286) && ($2 >= 12986) && ($2 <= 13286)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site2.matrix
awk '($1 >= 12986) && ($1 <= 13286) && ($2 >= 12986) && ($2 <= 13286)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site2.matrix

awk '($1 >= 36117) && ($1 <= 36417) && ($2 >= 36117) && ($2 <= 36417)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site3.matrix
awk '($1 >= 36117) && ($1 <= 36417) && ($2 >= 36117) && ($2 <= 36417)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site3.matrix
awk '($1 >= 36117) && ($1 <= 36417) && ($2 >= 36117) && ($2 <= 36417)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site3.matrix

awk '($1 >= 56733) && ($1 <= 57033) && ($2 >= 56733) && ($2 <= 57033)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site4.matrix
awk '($1 >= 56733) && ($1 <= 57033) && ($2 >= 56733) && ($2 <= 57033)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site4.matrix
awk '($1 >= 56733) && ($1 <= 57033) && ($2 >= 56733) && ($2 <= 57033)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site4.matrix

awk '($1 >= 58094) && ($1 <= 58394) && ($2 >= 58094) && ($2 <= 58394)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site5.matrix
awk '($1 >= 58094) && ($1 <= 58394) && ($2 >= 58094) && ($2 <= 58394)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site5.matrix
awk '($1 >= 58094) && ($1 <= 58394) && ($2 >= 58094) && ($2 <= 58394)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site5.matrix

awk '($1 >= 107479) && ($1 <= 107779) && ($2 >= 107479 ) && ($2 <= 107779)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site6.matrix
awk '($1 >= 107479) && ($1 <= 107779) && ($2 >= 107479 ) && ($2 <= 107779)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site6.matrix
awk '($1 >= 107479) && ($1 <= 107779) && ($2 >= 107479 ) && ($2 <= 107779)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site6.matrix

awk '($1 >= 111667) && ($1 <= 111967) && ($2 >= 111667) && ($2 <= 111967)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site7.matrix
awk '($1 >= 111667) && ($1 <= 111967) && ($2 >= 111667) && ($2 <= 111967)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site7.matrix
awk '($1 >= 111667) && ($1 <= 111967) && ($2 >= 111667) && ($2 <= 111967)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site7.matrix

awk '($1 >= 116234) && ($1 <= 116534) && ($2 >= 116234) && ($2 <= 116534)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site8.matrix
awk '($1 >= 116234) && ($1 <= 116534) && ($2 >= 116234) && ($2 <= 116534)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site8.matrix
awk '($1 >= 116234) && ($1 <= 116534) && ($2 >= 116234) && ($2 <= 116534)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site8.matrix

awk '($1 >= 116285) && ($1 <= 116585) && ($2 >= 116285) && ($2 <= 116585)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site9.matrix
awk '($1 >= 116285) && ($1 <= 116585) && ($2 >= 116285) && ($2 <= 116585)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site9.matrix
awk '($1 >= 116285) && ($1 <= 116585) && ($2 >= 116285) && ($2 <= 116585)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site9.matrix

awk '($1 >= 118656) && ($1 <= 118956) && ($2 >= 118656) && ($2 <= 118956)' hap1_clone21_hic_20000_iced.matrix > hap1_clone21_hic_20000_iced_site10.matrix
awk '($1 >= 118656) && ($1 <= 118956) && ($2 >= 118656) && ($2 <= 118956)' hap1_clone25_hic_20000_iced.matrix > hap1_clone25_hic_20000_iced_site10.matrix
awk '($1 >= 118656) && ($1 <= 118956) && ($2 >= 118656) && ($2 <= 118956)' hap1_wt_hic_20000_iced.matrix > hap1_wt_hic_20000_iced_site10.matrix
```

Generate Hi-C contact maps using lib5c: [paper](https://www.sciencedirect.com/science/article/pii/S2405471219300675) and [package](https://bitbucket.org/creminslab/lib5c/src/master/).

After installing lib5c, the following data sets/dependencies can be downloaded:

1. [Hi-C: insertion-proximal matrcies](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137372)
2. [CTCF/RAD21 ChIP-seq: insertion-proximal bedgraphs](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152721)
3. [RNA-seq: insertion-proximal bedgraphs](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137374)
4. [Compartment eigenvalues](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137372)
5. [CTCF motifs](https://upenn.box.com/s/2zijrporwjjzeg7cohfoleoj2v5exv15)
6. [hap1_wt_hic_20000_ord.bed(genomic coordinates to bins)](https://upenn.box.com/s/en9396j6vkkxuz7g0ai0t4ebdp7kmf1a)
7. conversion.py

Under Python 2.7.11+:

``` python
# import modules
from __future__ import division
import pandas as pd
import numpy as np
import scipy.sparse as sparse
import matplotlib as mpl
import matplotlib.colors as colors

from lib5c.contrib.pybigwig.bigwig import BigWig
from lib5c.util.mathematics import symmetrize
from lib5c.parsers.bed import load_features
from lib5c.contrib.interlap.util import features_to_interlaps, query_interlap

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from lib5c.plotters.extendable import ExtendableHeatmap
from conversion import grange_to_slice, slice_to_grange

import os
import glob

# update working directory
os.chdir('<working directory>')

# choose to plot for C21 insertions vs. C25 insertions:
clone_for_rna = "C21S"

# set heatmap ranges and upper/lower bounds for contacts
bin_each_side = [25, 30, 35, 40, 50, 60, 70, 80, 100]
max_cutoff = [15, 20, 25, 30, 35, 40]
min_cutoff = 0.1043
chip_size = '7%'


ctcf_motifs = load_features('pwmscan_hg19_MA01391_CTCF_5E-5.bed')

# assign corresponding data sets to each Hi-C sample
for file in glob.glob("*.matrix"):
	matrix_name = file
	sample = file.split("_")[1]
	if sample == "clone21":
		sf = 1
		rna_seq_id = "1734"
		ctcf_id = "2597"
		rad21_id = "2598"
	elif sample == "clone25":
		sf = 1.0472
		rna_seq_id = "1736"
		ctcf_id = "2601"
		rad21_id = "2602"
	elif sample == "wt":
		sf = 1.0505
		rna_seq_id = "1732"
		ctcf_id = "2593"
		rad21_id = "2594"
	elif sample == "sc79":
		sf = 160774393/161448929
		rna_seq_id = "2564"
		ctcf_id = "2605"
		rad21_id = "2606"
	elif sample == "c21-2-13":
		sf = 168920452/161448929
		rna_seq_id = "2694"
		ctcf_id = "2729"
		rad21_id = "2730"
	elif sample == "sc79-2-3":
		sf = 169715660/161448929
		rna_seq_id = "2697"
		ctcf_id = "2734"
		rad21_id = "2735"
	elif sample == "wt-2kb":
		sf = 169060033/161448929
		ctcf_id = "2724"
		rad21_id = "2725"

	# obtain info on which insertion site to plot
	site = (file.split("_")[-1]).split(".")[0]
	s = site[-1]
	if s == "0":
		s = "10"
	if s == "o":
		s = "endo"
		clone_for_rna= ""

	# match RNA-seq data
	if rna_seq_id < 1740:
		rna_seq_plus_track = load_features(rna_seq_id + ".hg19.plus.bw." + clone_for_rna + s + ".bedGraph")
		rna_seq_minus_track = load_features(rna_seq_id + ".hg19.minus.bw." + clone_for_rna + s + ".bedGraph")
	else:
		rna_seq_plus_track = load_features(rna_seq_id + ".plus.hg19.bw." + clone_for_rna + s + ".bedGraph")
		rna_seq_minus_track = load_features(rna_seq_id + ".minus.hg19.bw." + clone_for_rna + s + ".bedGraph")

	# match chip tracks
	ctcf_track = load_features(ctcf_id + ".hg19.cpm.bw." + clone_for_rna + s + ".bedGraph")
	rad21_track = load_features(rad21_id + ".hg19.cpm.bw." + clone_for_rna + s + ".bedGraph")
	chip_peaks = load_features(ctcf_id + ".hg19.sig20cutoff.broadpeak")

	# match compartment eigenvalues
	with open(ctcf_id + ".hg19.cpm.bw." + clone_for_rna + s + ".bedGraph") as f_for_chrom:
		first_line_for_chrom = f_for_chrom.readline()
	if first_line_for_chrom[4] == '\t':
		chr_for_eigen = first_line_for_chrom[3]
	else:
		chr_for_eigen = first_line_for_chrom[3:5]
	eigen_track = load_features("hap1_" + sample + "_eigen_" + chr_for_eigen + ".bedgraph")

	# plot for each distance range & max cutoff
	for bin in bin_each_side:
		for ma in max_cutoff:
			file_name = sample + '_' + site + '_' + str(2*bin) + '_' + 'bins' + '_' + 'cutoff' + str(ma) + '_' + chip_size + "." + "png"
			df = pd.read_csv(matrix_name, sep='\t', names=['i', 'j', 'v'])
			csr = sparse.coo_matrix((df.v, (df.i-1, df.j-1))).tocsr()
			start_bin = min(df.i) - 1
			end_bin = max(df.i)
			left_bin = start_bin + 150 - bin
			right_bin = end_bin - 150 + bin
			s = slice(left_bin, right_bin)
			grange = slice_to_grange(s)
			h = ExtendableHeatmap(
			    symmetrize(csr[s, s].toarray().T)/sf,
			    grange_x=grange,
			    colormap='red8',
			    colorscale=(min_cutoff, ma)
			)
			cbar = h.add_colorbar()
			cbar.set_label('Normalizated Interaction Count')
			peaks = features_to_interlaps(chip_peaks)[grange['chrom']]
			motifs = ctcf_motifs[grange['chrom']]
			h.add_refgene_stacks('hg19')
			h.add_chipseq_tracks(rna_seq_plus_track[h.grange_x['chrom']], size=chip_size, axis_limits=(0, 1), name='rnaseq_plus')
			h.add_chipseq_tracks(rna_seq_minus_track[h.grange_x['chrom']], size=chip_size, axis_limits=(-1, 0), name='rnaseq_minus')
			h.add_chipseq_tracks(eigen_track[h.grange_x['chrom']], size='10%', axis_limits=(-0.04, 0.04), name='eigen')
			h.add_motif_tracks([m for m in motifs if query_interlap(peaks, m)], colors={'+': 'r', '-': 'b'}, size=chip_size)
			h.add_chipseq_tracks(ctcf_track[h.grange_x['chrom']], size=chip_size, axis_limits=(0, 6), name='ctcf')
			h.add_chipseq_tracks(rad21_track[h.grange_x['chrom']], size=chip_size, axis_limits=(0, 3), name='rad21')
   			h['vertical_eigen'].set_xticks([])
   			h['horizontal_eigen'].set_yticks([])
   			h['vertical_rnaseq_plus'].set_xticks([])
   			h['horizontal_rnaseq_plus'].set_yticks([])
   			h['vertical_rnaseq_minus'].set_xticks([])
   			h['horizontal_rnaseq_minus'].set_yticks([])
   			h['vertical_ctcf'].set_xticks([])
   			h['horizontal_ctcf'].set_yticks([])
   			h['vertical_rad21'].set_xticks([])
   			h['horizontal_rad21'].set_yticks([])
			h.fig.savefig(file_name, dpi=800, bbox_inches='tight')
			mpl.pyplot.close('all')

```

### Log2 comparison maps
Similar as above, using lib5c under Python 2.7.11+:

``` python
# import modules
import pandas as pd
import numpy as np
import scipy.sparse as sparse
import matplotlib as mpl
import matplotlib.colors as colors

from lib5c.contrib.pybigwig.bigwig import BigWig
from lib5c.util.mathematics import symmetrize
from lib5c.parsers.bed import load_features
from lib5c.contrib.interlap.util import features_to_interlaps, query_interlap

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from lib5c.plotters.extendable import ExtendableHeatmap
from conversion import grange_to_slice, slice_to_grange

import os
import glob

# specify samples
insert_clone = "clone25"
sample_2 = "wt"
sample_3 = "clone21"

# specify distance ranges and log2 fold cutoff
bin_each_side = [20, 50, 80]
cutoff = [1.5]

# generate pairwise log2 comparison maps
for file in glob.glob("hap1*" + insert_clone + "*.matrix"):
	site = (file.split("_")[-1]).split(".")[0]
	matrix_1_name = file
	sample_1 = file.split("_")[0]
	df_1 = pd.read_csv(matrix_1_name, sep='\t', names=['i', 'j', 'v'])
	csr_1 = sparse.coo_matrix((df_1.v, (df_1.i-1, df_1.j-1))).tocsr()
	matrix_2_name = 'hap1' + '_' + sample_2 + '_hic_20000_iced_' + file.split("_")[-1]
	df_2 = pd.read_csv(matrix_2_name, sep='\t', names=['i', 'j', 'v'])
	csr_2 = sparse.coo_matrix((df_2.v, (df_2.i-1, df_2.j-1))).tocsr()
	matrix_3_name = 'hap1' + '_' + sample_3 + '_hic_20000_iced_' + file.split("_")[-1]
	df_3 = pd.read_csv(matrix_3_name, sep='\t', names=['i', 'j', 'v'])
	csr_3 = sparse.coo_matrix((df_3.v, (df_3.i-1, df_3.j-1))).tocsr()
	for bin in bin_each_side:
		start_bin = min(df_1.i) - 1
		end_bin = max(df_1.i)
		left_bin = start_bin + 150 - bin
		right_bin = end_bin - 150 + bin
		s = slice(left_bin, right_bin)
		grange = slice_to_grange(s)
		for co in cutoff: 
			file_name = insert_clone + '_' + site + '_' + insert_clone + '_' + 'vs' + '_' + sample_2 + '_' + str(2*bin) + '_' + 'bins' + '_' + str(co) + 'cutoff' + "." + "png"
			h = ExtendableHeatmap(
			    np.log2((symmetrize(csr_1[s, s].toarray().T)/1)/(symmetrize(csr_2[s, s].toarray().T)/1.0505)),
			    grange_x=grange,
			    colormap='RdBu_r',
			    colorscale=(-co, co)
			)
			cbar = h.add_colorbar()
			cbar.set_label('log2 fold change')
			h.fig.savefig(file_name, dpi=800, bbox_inches='tight')
			file_name = insert_clone + '_' + site + '_' + insert_clone + '_' + 'vs' + '_' + sample_3 + '_' + str(2*bin) + '_' + 'bins' + '_' + str(co) + 'cutoff' + "." + "png"
			h = ExtendableHeatmap(
				np.log2((symmetrize(csr_1[s, s].toarray().T)/1)/(symmetrize(csr_3[s, s].toarray().T)/1.0472)),
			    grange_x=grange,
			    colormap='RdBu_r',
			    colorscale=(-co, co)
			)
			cbar = h.add_colorbar()
			cbar.set_label('log2 fold change')
			h.fig.savefig(file_name, dpi=800, bbox_inches='tight')
			file_name = insert_clone + '_' + site + '_' + sample_3 + '_' + 'vs' + '_' + sample_2 + '_' + str(2*bin) + '_' + 'bins' + '_' + str(co) + 'cutoff' + "." + "png"
			h = ExtendableHeatmap(
				np.log2((symmetrize(csr_3[s, s].toarray().T)/1.0472)/(symmetrize(csr_2[s, s].toarray().T)/1.0505)),
			    grange_x=grange,
			    colormap='RdBu_r',
			    colorscale=(-co, co)
			)
			cbar = h.add_colorbar()
			cbar.set_label('log2 fold change')
			h.fig.savefig(file_name, dpi=800, bbox_inches='tight')
			mpl.pyplot.close('all')

```

### Plotting insulation scores for insertion-proximal regions
Requires:

1. sparse_to_dense.r
2. mat_to_insulation.r
3. plot_columns.r
4. insu_across_insertions.r

In R, take Clone 21 Site 4 (C21S4) for example. 

Load functions:
``` r
source("sparse_to_dense.R")
source("mat_to_insulation.R")
source("plot_columns.R")
source("insu_across_insertions.R")
```

Read extracted matrices and convert them from sparse to dense: 
``` r
C21_C21S4 <- read.table('hap1_clone21_hic_20000_iced_site4.matrix')
C21_forC21site4 <- spar.to.dense(C21_C21S4)

WT_C21S4 <- read.table('hap1_wt_hic_20000_iced_site4.matrix')
WT_forC21site4 <- spar.to.dense(WT_C21S4)

C25_C21S4 <- read.table('hap1_clone25_hic_20000_iced_site4.matrix')
C25_forC21site4 <- spar.to.dense(C25_C21S4)
```

Scale/normalize to account for minor coverage differences: division by (no. of valid pairs/approx. # of genomic bins).
``` r
C21_forC21site4_normed <- C21_forC21site4/(161448929/154770)
WT_forC21site4_normed <- WT_forC21site4/(169604428/154770)
C25_forC21site4_normed <- C25_forC21site4/(169069685/154770)
```

Compute insulation matrix:
``` r
C21_forC21site4_insu <- mat.to.insulation(C21_forC21site4_normed,30)
C25_forC21site4_insu <- mat.to.insulation(C25_forC21site4_normed,30)
WT_forC21site4_insu <- mat.to.insulation(WT_forC21site4_normed,30)
```

Generate insulation scores of the specified distance range, througout a series of window sizes.  
``` r
insulation_across_insertions(WT_forC21site4_insu,C25_forC21site4_insu,C21_forC21site4_insu,150,50,30,"C21_site4_IS_50bins.pdf")
```
