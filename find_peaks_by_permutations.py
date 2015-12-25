"""
Load txpt objects and read locations within txpts.

For each txpt: call peaks by read shuffling (in negative and positive) and by CWT.

Rank by p value for the shuffle.

Apply a BH FDR.

INPUT:
Config file.

OUTPUT:
Outputs to lib['clusters_dir'] from the config.ini file.

Run as;
python find_peaks_by_permutations.py -c <directory with config.py> --use lib

From config.ini, this uses:
lib[control_bed1]...
lib[expt_bed1]...
lib['gtf_raw']
lib['read_beds']
lib['gtf_with_names']
exp_bedgraph_minus
exp_bedgraph_plus

cluster_tables_to_peak_files:
gtf_with_names
gtf_raw
clusters_dir
permutation_peaks_dir

subpeaks_ui.py
lib['coverage_wigs.]
exp_bedgraph_minus
exp_bedgraph_plus
lib['chr_sizes']

gtf_pickle
feat_loc_dir
gtf_one_txpt_per_gene
exp_bedgraph_minus
exp_bedgraph_plus
coverage_wigs = bedgraphs_folder
txpt_obj_path
bedgraphs_folder
"""

import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import numpy as np
import pandas
import collections
import re

import cluster_reads
import p_values_of_clusters


def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument('-c', '--config_dir',
                        default='config/',
                        help='''Directory holding a config.py file.''')
    parser.add_argument('-v', '--no_ui',
                        help='No user input.',
                        default=False, action='store_true')
    args = parser.parse_args()
    return args


def get_bed_files(bed_folder='beds/', gtf=None, args=None, lib=None):
    ht_exons = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    gtfd = collections.defaultdict(set)
    for feature in gtf:
        if feature.type == 'exon':
            ht_exons[feature.iv] += feature.name
            gtfd[feature.name].add(feature)
    control_names = [lib[x] for x in lib.keys() if re.match('control_bed.*', x) is not None]
    exp_names = [lib[x] for x in lib.keys() if re.match('exp_bed.*', x) is not None]
    all_names = control_names + exp_names
    bed_file_list = [
        "/".join([bed_folder, x]) for x in control_names + exp_names
    ]
    reads_by_gene = {}
    counts_by_gene = {}
    ga = {}
    for fname in bed_file_list:
        (_ga, reads_by_gene[fname], counts_by_gene[fname]) = read_a_bed_file(fname, ht_exons)
        ga[fname] = _ga
    return ga, reads_by_gene, counts_by_gene


def read_a_bed_file(fname, ht_exons):
    import time
    start_time = time.time()
    counts = collections.defaultdict(int)
    reads_by_gene = collections.defaultdict(set)
    _ga = HTSeq.GenomicArray('auto', stranded=True)
    if not os.path.exists(fname):
        print "bed file not found: {o}".format(o=fname)
        return counts, reads_by_gene
    for n, line in enumerate(open(fname)):
        # if n > 1e6: break
        if not n % 1e6:
            min_per_line = float(time.time() - start_time)/float(60 * max([1., n]))
            print "Line {n}: Elapsed: {te} s. Minutes for 10 million reads: {f}".format(
                n=n, te=float(time.time() - start_time),
                f=min_per_line * 1e7/60.
            )
        s = line.rstrip('\n').split('\t')
        read_iv = HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5])
        _ga[read_iv] += 1
        gene_ids = set()
        for iv, val in ht_exons[read_iv].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[gene_id] += 1
            reads_by_gene[gene_id].add(read_iv)
    print "Stats on {a}:".format(a=fname)
    cluster_reads.stats(reads_by_gene)
    return _ga, reads_by_gene, counts


def dir_creation(lib):
    dir_list = ['bedgraphs/', 'bed_uncollapsed/' 'bed_collapsed/',
                'figs/', 'logs/']
    for adir in dir_list:
        if not os.path.exists(lib['top'] + '/' + adir):
            os.system('mkdir ' + lib['top'] + '/' + adir)
        if not os.path.exists(lib['clusters_dir']):
            os.system('mkdir ' + lib['clusters_dir'])


if __name__ == '__main__':
    args = read_args()
    print args
    sys.path.insert(0, args.config_dir)
    import config
    lib = config.config()
    del sys.path[0]
    dir_creation(lib)
    logging.basicConfig(
        filename='logs/%s_find_peaks_by_perm.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    # Load raw reads.
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    reads_by_gene, counts_by_gene = get_bed_files(
        bed_folder=lib['read_beds'], gtf=gtf,
        args=args, lib=lib)
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf_with_names'])
    if args.no_ui:
        cluster_reads.run(
            reads_by_gene, counts_by_gene, exons_as_rows,
            lib=lib, args=args)
        print "Successfully finished."
        sys.exit()
    while True:
        try:
            cluster_reads.run(
                reads_by_gene, counts_by_gene, exons_as_rows,
                lib=lib, args=args)
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(config)
                reload(cluster_reads)
                reload(p_values_of_clusters)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

