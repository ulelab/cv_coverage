"""Comparative visualisation of specific motifs across data or regions

The user should define:
1. List of 1 or more kmers that will be visualized as a group (based on average counts)
2. List of data (1 or more) to compare. Data is in form of bed files containing crosslinks.
3. List of  genomic regions (1 or more) to compare. possible regions are 'genome', 'intergenic',
'intron', 'UTR3', 'other_exon'(comprised of 'UTR5' and 'CDS') and 'ncRNA'.

First step is regional thresholding to obtain thresholded crosslinks. This approach takes crosslinks
in all peaks within a region to define threshold and so introduces an element of intra-regional comparison.
Regions for thresholding as defined in the following way:
- all exons in the same gene (5’UTR, CDS, 3’UTR, or all exons in ncRNAs) are considered one region,
- each intron is its own region,
- each intergenic region is its own region.

Draw the selected kmer positional distribution around thresholded xcrosslinks (-window..100) for the
multiple data/regions as defined by a user. If multiple data AND regions are given,
then each region is shown on a separate plot.
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools as pbt
from random import randint
from itertools import islice
from multiprocessing import  Pool
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


REGION_SITES = {
    'genome': ['intron', 'CDS', 'UTR3', 'UTR5', 'ncRNA', 'intergenic'],
    'whole_gene': ['intron', 'CDS', 'UTR3', 'UTR5'],
    'intergenic': ['intergenic'],
    'intron': ['intron'],
    'ncRNA': ['ncRNA'],
    'other_exon': ['UTR5', 'CDS'],
    'UTR3': ['UTR3'],
    'UTR5': ['UTR5']
}
REGIONS_QUANTILE = [
    'intron',
    'intergenic',
    'cds_utr_ncrna',
]
REGIONS_MAP = {}
TEMP_PATH = None

# overriding pybedtools to_dataframe method to avoid warning
def to_dataframe_fixed(self, *args, **kwargs):
    """
    Create a pandas.DataFrame, passing args and kwargs to pandas.read_csv.

    This function overrides pybedtools function to avoid FutureWarning:
    read_table is deprecated, use read_csv instead... Pandas must be
    imported as pd, it is advisable to specify dtype and names as well.
    """
    return pd.read_csv(self.fn, header=None, sep='\t', *args, **kwargs)


pbt.BedTool.to_dataframe = to_dataframe_fixed  # required for overriding


def get_name(s_file):
    """Return sample name from file path."""
    return s_file.split('/')[-1].replace('.gz', '').replace('.bed', "").replace('.xl', "")



def parse_region_to_df(region_file):
    """Parse GTF to pandas.DataFrame."""
    return pd.read_csv(
        region_file,
        names=['chrom', 'second', 'region', 'start', 'end', 'sixth', 'strand', 'eighth', 'id_name_biotype'],
        sep='\t',
        header=None,
        dtype={
            'chrom': str, 'second': str, 'region': str, 'start': int, 'end': int, 'sixth': str, 'strand': str,
            'eight': str, 'id_name_biotype': str
        }
    )


def filter_cds_utr_ncrna(df_in):
    """Filter regions CDS, UTR5, UTR3 and ncRNA by size and trim."""
    utr5 = df_in.region == 'UTR5'
    cds = df_in.region == 'CDS'
    utr3 = df_in.region == 'UTR3'
    ncrna = df_in.region == 'ncRNA'
    size = df_in.end - df_in.start >= 100
    df_out = df_in[(utr5 & size) | (cds & size) | (utr3 & size) | ncrna].copy()
    df_out.loc[df_out['region'] == 'CDS', ['start']] = df_out.start + 30
    df_out.loc[df_out['region'] == 'CDS', ['end']] = df_out.end - 30
    return df_out


def filter_intron(df_in, min_size):
    """Filter intron regions to remove those smaller than min_size."""
    # remove regions shorter then min_size
    df_out = df_in[df_in.end - df_in.start >= min_size].copy()
    return df_out


def get_regions_map(regions_file):
    """Prepare temporary files based on GTF file that defines regions."""
    df_regions = pd.read_csv(
        regions_file, sep='\t', header=None,
        names=['chrom', 'second', 'region', 'start', 'end', 'sixth', 'strand', 'eighth', 'id_name_biotype'],
        dtype={
            'chrom': str, 'second': str, 'region': str, 'start': int, 'end': int, 'sixth': str, 'strand': str,
            'eight': str, 'id_name_biotype': str})
    df_intergenic = df_regions.loc[df_regions['region'] == 'intergenic']
    df_cds_utr_ncrna = df_regions.loc[df_regions['region'].isin(['CDS', 'UTR3', 'UTR5', 'ncRNA'])]
    df_intron = df_regions.loc[df_regions['region'] == 'intron']
    df_cds_utr_ncrna = filter_cds_utr_ncrna(df_cds_utr_ncrna)
    df_intron = filter_intron(df_intron, 100)
    to_csv_kwrgs = {'sep': '\t', 'header': None, 'index': None}
    df_intron.to_csv('{}/intron_regions.bed'.format(TEMP_PATH), **to_csv_kwrgs)
    df_intergenic.to_csv('{}/intergenic_regions.bed'.format(TEMP_PATH), **to_csv_kwrgs)
    df_cds_utr_ncrna.to_csv('{}/cds_utr_ncrna_regions.bed'.format(TEMP_PATH), **to_csv_kwrgs)

def get_sequences(sites, fasta, fai, window_l, window_r, merge_overlaps=False):
    """Get genome sequences around positions defined in sites."""
    sites = pbt.BedTool(sites).sort()
    sites_extended = sites.slop(l=window_l, r=window_r, g=fai)  # noqa
    if merge_overlaps:
        sites_extended = sites_extended.merge(s=True)
    seq_tab = sites_extended.sequence(s=True, fi=fasta, tab=True)
    return [line.split("\t")[1].strip() for line in open(seq_tab.seqfn)]


def intersect(interval_file, s_file):
    """Intersect two BED files and return resulting BED file."""
    if interval_file:
        result = pbt.BedTool(s_file).intersect(
            pbt.BedTool(interval_file), s=True,
            nonamecheck=True,
        ).saveas()
    else:
        result = pbt.BedTool(s_file)
    if len(result) >= 1:
        return result


def intersect_merge_info(region, s_file):
    """Intersect while keeping information from region file."""
    interval_file = REGIONS_MAP[region]
    try:
        df_1 = intersect(interval_file, s_file).to_dataframe(
            names = ['chrom', 'start', 'end', 'name', 'score', 'strand'],
            dtype={'chrom': str, 'start': int, 'end': int, 'name': str, 'score': float, 'strand': str}
        )
        df_1 = df_1.groupby(['chrom', 'start', 'end', 'strand'], as_index=False)['score'].sum(axis=0)
        df_1['name'] = '.'
        df_2 = intersect(s_file, interval_file).to_dataframe(
            names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'],
            dtype={'seqname': str, 'source': str, 'feature': str, 'start': int, 'end': int, 'score': str,
                'strand': str, 'frame': str, 'attributes': str}
        )
        df_2.drop_duplicates(subset=['seqname', 'start', 'end', 'strand'], keep='first')
    except AttributeError:
        print(f'{s_file} might have no sites in {region}')
        return
    df_2 = df_2.drop(columns=['source', 'score', 'frame', 'start']).rename(index=str, columns={"seqname": "chrom"})
    return pd.merge(df_1, df_2, on=['chrom', 'strand', 'end'])

def cut_per_chrom(chrom, df_p, df_m, df_peaks_p, df_peaks_m):
    """Split data by strand then apply pandas cut to each strand.

    Pandas cut uses IntervalIndex (done from the peaks file) to
    assign each site its peak. Finally merges strands.
    """
    df_temp_p = df_peaks_p[df_peaks_p['chrom'] == chrom].copy()
    df_temp_m = df_peaks_m[df_peaks_m['chrom'] == chrom].copy()
    df_xl_p = df_p[df_p['chrom'] == chrom].copy()
    df_xl_m = df_m[df_m['chrom'] == chrom].copy()
    left_p = np.array(df_temp_p['start'])
    right_p = np.array(df_temp_p['end'])
    left_m = np.array(df_temp_m['start'])
    right_m = np.array(df_temp_m['end'])
    interval_index_p = pd.IntervalIndex.from_arrays(left_p, right_p, closed='left')
    interval_index_m = pd.IntervalIndex.from_arrays(left_m, right_m, closed='left')
    df_xl_p['cut'] = pd.cut(df_xl_p['start'], interval_index_p)
    df_xl_m['cut'] = pd.cut(df_xl_m['start'], interval_index_m)
    return pd.concat([df_xl_p, df_xl_m], ignore_index=True)


def cut_sites_with_region(df_sites, df_region):
    """Find peak interval the crosslinks belong to."""
    df_p = df_sites[df_sites['strand'] == '+'].copy()
    df_m = df_sites[df_sites['strand'] == '-'].copy()
    df_region_p = df_region[df_region['strand'] == '+'].copy()
    df_region_m = df_region[df_region['strand'] == '-'].copy()
    df_cut = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes', 'cut'])
    for chrom in set(df_region['chrom'].values):
        df_temp = cut_per_chrom(chrom, df_p, df_m, df_region_p, df_region_m)
        df_temp = df_temp[df_cut.columns]
        df_cut = pd.concat([df_cut, df_temp], ignore_index=True)
    return df_cut.dropna(axis=0)


def percentile_filter_xlinks(df_in, percentile=0.7):
    """Calculate threshold and filter sites by it."""
    df_in['cut'] = df_in['cut'].astype(str)
    df_in['quantile'] = df_in['cut'].map(df_in.groupby('cut').quantile(q=percentile)['score'])
    df_in = df_in[df_in['score'] > df_in['quantile']]
    return df_in[['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes']]


def get_threshold_sites(s_file, percentile=0.7):
    """Apply crosslink filtering based on dynamical thresholds.

    Regions for thresholds are defined as follows: introns and
    intergenic regions are each idts own region, for CDS, UTR and ncRNA
    each gene is a region. After region determination threshold based on
    percentile are applied and finally threshold crosslinks sites are
    sorted.
    """
    df_out = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes'])
    for region in REGIONS_QUANTILE:
        print(f'Thresholding {region}')
        df_reg = intersect_merge_info(region, s_file)
        print(f'lenght of df_reg for {region} is: {len(df_reg)}')
        if df_reg is None:
            continue

        if region == 'cds_utr_ncrna':
            df_reg.name = df_reg.attributes.map(lambda x: x.split(';')[1].split(' ')[1].strip('"'))
            df_reg['quantile'] = df_reg['name'].map(df_reg.groupby(['name']).quantile(q=percentile)['score'])
            df_filtered = df_reg[df_reg['score'] > df_reg['quantile']].drop(columns=['quantile'])
            df_out = pd.concat([df_out, df_filtered], ignore_index=True, sort=False)
        if region in ['intron', 'intergenic']:
            df_region = parse_region_to_df(REGIONS_MAP[region])
            df_cut = cut_sites_with_region(df_reg, df_region)
            df_filtered = percentile_filter_xlinks(df_cut)
            df_out = pd.concat([df_out, df_filtered], ignore_index=True, sort=False)
    return df_out.sort_values(by=['chrom', 'start', 'strand'], ascending=[True, True, True]).reset_index(drop=True)




# def get_all_sites(s_file):
#     """Get crosslink data into appropriate dataframe without thresholding."""
#     df_out = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes'])
#     for region in REGIONS_QUANTILE:
#         df_reg = intersect_merge_info(region, s_file)
#         if df_reg is None:
#             return
#         if df_reg.empty:
#             continue
#         if region == 'cds_utr_ncrna':
#             df_reg.name = df_reg.attributes.map(lambda x: x.split(';')[1].split(' ')[1].strip('"'))
#             df_reg['quantile'] = None
#             df_out = pd.concat([df_out, df_reg], ignore_index=True, sort=False)
#         if region in ['intron', 'intergenic']:
#             df_region = parse_region_to_df(REGIONS_MAP[region])
#             df_cut = cut_sites_with_region(df_reg, df_region)
#             df_filtered = df_cut[['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes']]
#             df_out = pd.concat([df_out, df_filtered], ignore_index=True, sort=False)
#     return df_out.sort_values(by=['chrom', 'start', 'strand'], ascending=[True, True, True]).reset_index(drop=True)


def get_all_sites(s_file):
    """Get crosslink data into appropriate dataframe without thresholding."""
    df_out = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes'])
    for region in REGIONS_QUANTILE:
        df_reg = intersect_merge_info(region, s_file)
        if df_reg is None:
            continue
        if df_reg.empty:
            continue
        if region == 'cds_utr_ncrna':
            df_reg.name = df_reg.attributes.map(lambda x: x.split(';')[1].split(' ')[1].strip('"'))
            df_reg['quantile'] = None
            df_out = pd.concat([df_out, df_reg], ignore_index=True, sort=False)
        if region in ['intron', 'intergenic']:
            df_region = parse_region_to_df(REGIONS_MAP[region])
            df_cut = cut_sites_with_region(df_reg, df_region)
            df_filtered = df_cut[['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature', 'attributes']]
            df_out = pd.concat([df_out, df_filtered], ignore_index=True, sort=False)
    return df_out.sort_values(by=['chrom', 'start', 'strand'], ascending=[True, True, True]).reset_index(drop=True)


def parallelize(func, sequences, scores, my_int1, my_int2, list_of_strings, use_scores, n_cores):
    split_seqs = [sequences[i::n_cores] for i in range(n_cores)]
    split_scores = [scores[i::n_cores] for i in range(n_cores)]
    pool = Pool(n_cores)
    iterable_args = [(split_seq, split_scores[i], my_int1, my_int2, list_of_strings, use_scores) for i, split_seq in enumerate(split_seqs)]
    results = pool.starmap(func, iterable_args)
    pool.close()
    pool.join()
    return results


def get_positions(s, kmer):
    def find_all(a_str, sub):
        """find indices of the substring in a string"""
        start = 0
        while True:
            start = a_str.find(sub, start)
            if start == -1: return
            yield start
            start += 1
    # for each sliding window we find which positions fit into each motif
    indices_extended = []
    indices = list(find_all(s, kmer))
    for i in indices:
        indices_extended.extend(range(i, (i + len(kmer))))
    position = [1 if x in set(indices_extended) else 0 for x in range(len(s))]
    # number of positions is summed up into score for the current window
    # score of the window with max score is returned, called h_max
    return position


def get_coverage_old(sequences, scores, window, kmer_length, kmer_list_input, use_scores):
    seq_pos_all_kmer_d = {k1: {k2: 0 for k2 in list(range(-(window + kmer_length), window + kmer_length + 1))} for k1 in kmer_list_input}
    for kmer in kmer_list_input:
        if not use_scores:
            seq_pos = [p for p in [get_positions(s, kmer) for s in sequences]]
            seq_pos = [l for l in seq_pos if len(l) and sum(l) > 0]
            seq_pos_sum = [sum(i) for i in zip(*seq_pos)]
        else:
            seq_pos = [[scores[i] * p for p in get_positions(s, kmer)] for i, s in enumerate(sequences)]
            seq_pos = [l for l in seq_pos if len(l) and sum(l) > 0]
            seq_pos_sum = [sum(i) for i in zip(*seq_pos)]
        for j, p in enumerate(seq_pos_sum):
            pos = j - (window + kmer_length)
            if p > 0:
                seq_pos_all_kmer_d[kmer][pos] += p
    seq_pos_all_kmer = {k: [u for u in v.values()] for k, v in seq_pos_all_kmer_d.items()}
    return seq_pos_all_kmer


def get_coverage(sequences, scores, window, kmer_length, kmer_list_input, use_scores):
    seq_pos_all = {k: 0 for k in list(range(-(window + kmer_length), window + kmer_length + 1))}
    for j, s in enumerate(sequences):
        seq_pos = [p for p in [get_positions(s, kmer) for kmer in kmer_list_input]]
        if not use_scores:
            seq_pos = [l for l in seq_pos if len(l) and sum(l) > 0]
            seq_pos_sum = [int(bool(sum(i))) for i in zip(*seq_pos)]
        else:
            seq_pos_sum = [int(bool(sum(i))) * scores[j] for i in zip(*seq_pos)]
        for j, p in enumerate(seq_pos_sum):
            pos = j - (window + kmer_length)
            if p > 0:
                seq_pos_all[pos] += p
    return seq_pos_all

def combine_results(results):
    combined_results = {k: 0 for k in results[0].keys()}
    for result in results:
        for k, v in result.items():
            combined_results[k] += v
    return combined_results


def combine_results_old(results):
    combined_results = {k: [] for k in results[0].keys()}
    for result in results:
        for k, v in result.items():
            combined_results[k].append(v)
    return {k: [sum(i) for i in zip(*v)] for k, v in combined_results.items()}


def run(sites_files_paths_list, kmer_list_input, region_list_input, kmer_length, genome, genome_fai,
regions_file, smoot, percentile=None, window=150, use_scores=False, n_cores=4, chunk_size=100000, cap=0):
    """Run the script"""
    global TEMP_PATH
    TEMP_PATH = os.environ['TMPDIR']
    os.makedirs('./results/', exist_ok=True)
    kmer_list_input = [kmer.replace('U', 'T') for kmer in kmer_list_input]
    txn_dict = {x: {y: 0 for y in region_list_input} for x in [get_name(s) for s in sites_files_paths_list]}
    df_out = pd.DataFrame()
    for sites_file in sites_files_paths_list:
        print(f'Analysing: {sites_file}')
        sample_name = get_name(sites_file)
        get_regions_map(regions_file)
        global REGIONS_MAP
        REGIONS_MAP = {
            'intron': '{}/intron_regions.bed'.format(TEMP_PATH),
            'intergenic': '{}/intergenic_regions.bed'.format(TEMP_PATH),
            'cds_utr_ncrna': '{}/cds_utr_ncrna_regions.bed'.format(TEMP_PATH)
        }
        df_params = pd.Series(data={"file_list": sites_files_paths_list,
                                "kmer_list": kmer_list_input,
                                "regions": region_list_input,
                                "kmer_length": kmer_length,
                                "genome": genome,
                                "genome_index_file": genome_fai,
                                "gtf_segmentation_file": regions_file,
                                "smoothing": smoot,
                                "percentile": percentile,
                                "window": window,
                                "use_scores": use_scores,
                                "n_cores": n_cores,
                                "chunk_size": chunk_size,
                                "cap" : cap,
                                })
        df_params.to_csv(f'./results/run_parameters.tsv', sep='\t', header=False)

        if percentile:
            print('Getting thresholded crosslinks')
            df_txn = get_threshold_sites(sites_file, percentile=percentile)
            print(f'{len(df_txn)} thresholded crosslinks')
        else:
            print('Getting all crosslinks')
            df_txn = get_all_sites(sites_file)
            if cap:
                df_txn['score'] = df_txn['score'].apply(lambda x: x if x <= cap else cap)


        for region in region_list_input:
            # Parse sites file and keep only parts that intersect with given region
            # for organisms which don't have introns, this will allow whole_gene analysis.
            regs_list = [r for r in REGION_SITES[region] if r in df_txn['feature'].unique().tolist()]
            if len(regs_list) >= 1:
                print(region, 'is represented by these annotated features:', regs_list)
            else:
                print(f'features corresponding to {region} were not found in regions gtf file.')
                # if this is the case remove this region from input
                region_list_input = [r for r in region_list_input if r!=region]
                continue
            df_sites = df_txn.loc[df_txn['feature'].isin(regs_list)]
            if percentile:
                print(f'{len(df_sites)} thresholded sites on {region}')
            else:
                print(f'{len(df_sites)} total sites on {region}')
            txn_dict[sample_name][region] = len(df_sites)
            scores = df_sites.score.values.astype(float)
            n = len(df_sites) // chunk_size + 1
            sites_chunks = np.array_split(df_sites, n)
            scores_chunk = np.array_split(scores, n)
            results_chunks = []
            for i, chunk in enumerate(sites_chunks):
                sites = pbt.BedTool.from_dataframe(chunk[['chrom', 'start', 'end', 'name', 'score', 'strand']])
                sequences = get_sequences(sites, genome, genome_fai, window + kmer_length, window + kmer_length)
                assert len(sequences) == len(scores_chunk[i])
                results = parallelize(get_coverage, sequences, scores_chunk[i], window, kmer_length, kmer_list_input, use_scores, n_cores)
                results_chunks.append(combine_results(results))
            combined_results = combine_results(results_chunks)
            if use_scores:
                total = scores.sum()
            else:
                total = len(df_sites)
            print('total ', total)
            seq_pos_all_kmer = {k: [v / total] for k, v in combined_results.items()}
            test_df = pd.DataFrame.from_dict(seq_pos_all_kmer).T
            index = pd.Index(test_df.index)
            # define name of column
            mean_col = '{}_{}'.format(sample_name, region)
            # calculate average of kmers of interest for each position and save results in new column
            # copy the column of average counts to dataframe used for plotting
            df_out[mean_col] = test_df.T.values[0]
            pbt.cleanup()
            print(f'{region} Finished')
    df_smooth = df_out.rolling(smoot, center=True, win_type='triang').mean()
    # slicing drops edge values that get NaN due to rolling mean
    df_smooth = df_smooth * 100
    print('All selected regions finished')
    # different cases of plotting the data, depending on how many plots need to be generated
    plot_list = []
    # if there is just one file we plot distributions for each region on the same plot
    if len(sites_files_paths_list) == 1:
        df_plot = df_smooth
        plot_list.append(df_plot)
    # if there are multiple file but single region we plot all files on same plot
    elif (len(sites_files_paths_list) > 1) and (len(region_list_input) == 1):
        df_plot = df_smooth
        plot_list.append(df_plot)
    # if there is more then one file and region we plot different files together for each region
    elif (len(sites_files_paths_list) > 1) and (len(region_list_input) > 1):
        for r in region_list_input:
            columns = [x for  x in df_smooth.columns if r in x]
            df_plot = df_smooth[columns]
            plot_list.append(df_plot)

    # save a tsv file used for plotting, which contains average counts of kmers of interest for each position
    kmer_list_name = [kmer.replace('T', 'U') for kmer in kmer_list_input]
    outfile_name = "_".join(kmer_list_name[:4])
    df_out.set_index(index, inplace=True)
    for region in region_list_input:
        df_out[[x for x in df_out.columns if region in x]].to_csv(
            f'./results/{sample_name}_{region}_{cap}_{outfile_name}.tsv',
            sep='\t',
            float_format='%.8f'
        )
    pd.DataFrame.from_dict(txn_dict, orient='index').to_csv(f'./results/no_of_txn_{sample_name}', sep='\t')
    print('Results saved')
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    lineplot_kwrgs = {'palette': "tab20", 'linewidth': 1, 'dashes': False, }
    for p in plot_list:
        p = p.set_index(index)
        plt.figure()
        sns_plot = sns.lineplot(data = p.loc[-window + kmer_length: window - kmer_length, :], **lineplot_kwrgs)
        plt.ylim(bottom=0)
        plt.ylabel('Coverage (%)')
        plt.xlabel('Positions of kmer start relative to crosslinks')
        legend = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
        if (len(region_list_input) == 1):
            region = p.columns[0].split('_')[-1]
            plt.title('Coverage of {} motif group in {}'.format(outfile_name, region))
            sns_plot.figure.savefig(
                f'./results/{sample_name}_positional_distribution_cap{cap}_{outfile_name}.pdf',
                bbox_extra_artists=(legend,),
                bbox_inches='tight'
            )
        else:
            region = p.columns[0].split('_')[-1]
            plt.title('Coverage of {} motif group in {}'.format(outfile_name, region))
            sns_plot.figure.savefig(
                f'./results/{sample_name}_{region}_positional_distribution_cap{cap}_{outfile_name}.pdf',
                bbox_extra_artists=(legend,),
                bbox_inches='tight'
            )
        plt.close()
    print('Plots saved')


if __name__ == "__main__":
    import sys

    xl_in = sys.argv[1].split(',')
    motifs = sys.argv[2].split(',')
    regions = sys.argv[3].split(',')
    kmer_len = int(sys.argv[4])
    fasta = sys.argv[5]
    fai = sys.argv[6]
    regions_file = sys.argv[7]
    smoothing = int(sys.argv[8])
    percentile = 'None' if 'None' else int(sys.argv[9])
    window = int(sys.argv[10])
    use_scores = sys.argv[11]
    n_cores = int(sys.argv[12])
    chunk_size = int(sys.argv[13])
    cap = int(sys.argv[14])

    if percentile.strip() == 'None':
        percentile = None
    else:
        percentile = int(percentile)

    if use_scores.strip() == 'True':
        use_scores = True
    else:
        use_scores = False

    run(xl_in, motifs, regions, kmer_len, fasta, fai, regions_file, smoothing, percentile,
        window, use_scores, n_cores, chunk_size, cap)








