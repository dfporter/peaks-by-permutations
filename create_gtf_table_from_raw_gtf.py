import re
import argparse
import sys


def add_wb_name_to_gtf(lib):
    fname = lib['gtf_raw']
    gtf = {}
    header = "0\t1\t2\t3\t4\t5\t6\t7\t8\tgene_name"
    header += "\ttranscript_id\ttranscript_name\texon_number\tgene_id\tbiotype\n"
    with open(fname, 'r') as f:
        for line in f:
            gene_id = re.match('.*gene_id "([^"]+)";.*', line)
            if gene_id is None: print line
            gene_id = gene_id.group(1).partition('.')[0]
            gene_name = re.match('.*gene_name "([^"]+)";.*', line)
            biotype = re.match('.*biotype "([^"]+)";.*', line)
            biotype = biotype.group(1)
            if gene_name in gtf:
                if gene_id == gtf[gene_name]['gene_id']:
                    gtf[gene_name]['lines'].append(line)
                else:
                    print "Two gene IDs for the same gene name:",
                    print "%s: %s vs %s" % (gene_name, gene_id, gtf[gene_name]['gene_id'])
            else:
                gtf[gene_name] = {'gene_id': gene_id, 'lines': [line], 'biotype': biotype}
    with open(lib['gtf'], 'w') as f:
        f.write(header)
        for gene in gtf:
            for line in gtf[gene]['lines']:
                f.write(line.rstrip('\n') + '\t' + '\t'.join([gtf[gene]['gene_id'], gtf[gene]['biotype']]) +'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Converts a raw gtf file into a table with \
        a header and gene_id and biotype columns.''')
    parser.add_argument(
        '-c', '--config_dir', default='analysis/',
        help='Directory of config.py file.'
    )
    args = parser.parse_args()
    sys.path.insert(0, args.config_dir)
    import config
    lib = config.config()
    add_wb_name_to_gtf(lib)