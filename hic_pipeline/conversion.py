import numpy as np
import pandas as pd


bin_df = pd.read_csv('hap1_wt_hic_20000_ord.bed', sep='\t',
                     names=['chrom', 'start', 'end', 'bin_idx'])
res = bin_df.loc[0, 'end'] - bin_df.loc[0, 'start']
chroms = bin_df['chrom'].unique()
offsets = bin_df.groupby('chrom', sort=False)['bin_idx'].min().values - 1
offsets_by_name = dict(zip(chroms, offsets))


def coord_to_bin(chrom, coord):
    return coord / res + offsets_by_name[chrom]


def bin_to_coord(bin_idx):
    chrom = chroms[np.argmin(bin_idx >= offsets) - 1]
    return chrom, (bin_idx - offsets_by_name[chrom]) * res


def grange_to_slice(grange):
    return slice(coord_to_bin(grange['chrom'], grange['start']),
                 coord_to_bin(grange['chrom'], grange['end']))


def slice_to_grange(bin_slice):
    chrom, start = bin_to_coord(bin_slice.start)
    _, end = bin_to_coord(bin_slice.stop)
    return {'chrom': chrom, 'start': start, 'end': end}
