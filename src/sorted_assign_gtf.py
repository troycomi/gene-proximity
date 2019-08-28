from argparse import ArgumentParser
import random


def main():
    random.seed(0)
    options = get_options()
    priority_order = {
        "intergenic": 1,
        "downstream": 2,
        "upstream": 3,
        "intron": 4,
        "3UTR": 5,
        "5UTR": 6,
        "exon": 7,
        "promoter": 8,  # highest priority
    }

    with open(options.gtf, 'r') as gtf, \
            open(options.locations, 'r') as locations, \
            open(options.output, 'w') as output, \
            open(options.priority, 'w') as priority:

        gtf_line = gtf.readline()
        loc_line = locations.readline()
        current_entry = None
        while gtf_line and loc_line:
            # remove any scaffolds that are unmapped
            if len(gtf_line.split()[0]) > 3:
                gtf_line = gtf.readline()
                continue

            split_gtf = gtf_line.split('\t')
            gtf_id = split_gtf[-1].split(';')[0].split('=')[1].strip('";')
            g_chrom = split_gtf[0]
            g_start = int(split_gtf[3])
            g_end = int(split_gtf[4])
            g_type = split_gtf[2]

            chrom, pos = loc_line.split()
            pos = int(pos)

            if g_chrom < chrom or pos > g_end:
                gtf_line = gtf.readline()
                continue

            if g_chrom > chrom or pos < g_start:
                loc_line = locations.readline()
                if current_entry:
                    priority.write(
                        f'{chrom}\t{pos}\t'
                        f'{random.sample(current_entry[1],1)[0]}\t'
                        f'{current_entry[0]}\n')
                continue

            # match!
            output.write(f'{chrom}\t{pos}\t{gtf_id}\t{g_type}\n')
            if not current_entry or \
                    priority_order.get(current_entry[0], 0) < \
                    priority_order.get(g_type, 0):
                current_entry = [g_type, {gtf_id}]
            else:
                current_entry[1].add(gtf_id)

            gtf_line = gtf.readline()


def get_options(arguments=None):
    '''read command line arguments and return the corresponding arguments'''
    parser = ArgumentParser(description='Assign ATACseq peaks to genes')

    parser.add_argument("-g", "--gtf",
                        help='sorted gtf file')
    parser.add_argument("-l", "--locations",
                        help='sorted list of chromosome/positions to find')
    parser.add_argument("-o", "--output",
                        help='all outputs')
    parser.add_argument("-p", "--priority",
                        help='just the highest priority feature')
    return parser.parse_args(arguments)


if __name__ == '__main__':
    main()
