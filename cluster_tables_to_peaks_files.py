"""
Load the clusters output by find_peaks_by_permutations.

"""
import HTSeq
import argparse
import os
import traceback
import sys
import pandas
import glob

import cluster_reads
import p_values_of_clusters
import cluster_combine

def read_args():
    parser = argparse.ArgumentParser(description='''
        Convert cluster tables to peaks files by comparing replicates/p values.''')
    parser.add_argument('-c', '--config_dir',
                        default='config/',
                        help='''Directory holding a config.py file.''')
    parser.add_argument('-v', '--no_ui',
                        help='No user input.',
                        default=False, action='store_true')

    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = read_args()
    print args
    sys.path.insert(0, args.config_dir)
    import config
    lib = config.config()
    del sys.path[0]
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    # get_txpts() uses the peaks argument to select what txpts to load by a call to
    #  set(peaks['gene_name'].tolist()). That is, peaks is a dataframe.
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf'])
    if args.no_ui:
        cluster_combine.run(args, lib, gtf, exons_as_rows)
        sys.exit()
    while True:
        try:
            cluster_combine.run(args, lib, gtf, exons_as_rows)
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
                reload(cluster_combine)
                reload(cluster_reads)
                reload(p_values_of_clusters)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

