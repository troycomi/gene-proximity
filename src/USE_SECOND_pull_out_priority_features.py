# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:31:05 2018

@author: Beryl Jones
"""

import random
from optparse import OptionParser
random.seed(a=1234)

parser = OptionParser()

parser.add_option('-i', '--input_file', dest='infile')
parser.add_option('-o', '--output_file', dest='outfile')

(options, args) = parser.parse_args()

# read in file
infile = options.infile
outfile = options.outfile

# chr	coord  proxtype	gene

priority_order = {
                   "intergenic": 1,
                   "downstream": 2,
                   "upstream": 3,
                   "intron": 4,
                   "three_utr": 5,
                   "five_utr": 6,
                   "exon": 7,
                   "promoter": 8,  # highest priority
                  }

with open(infile, 'r') as reader:
    with open(outfile, 'w') as writer:
        writer.write(reader.readline())

        current_peak = ""
        highest_priority = ""
        gene_list = []
        for l in reader.readlines():
            tokens = l.strip().split('\t')
            peak = "\t".join(tokens[0:2])
            if peak != current_peak:
                if current_peak == "":
                    pass  # do nothing, first peak
                else:
                    writer.write("{}\t{}\t{}\n".format(
                        current_peak,
                        highest_priority,
                        random.choice(gene_list)))
                # first entry
                current_peak = peak
                highest_priority = tokens[2]
                gene_list = [tokens[3]]
            else:
                if priority_order[highest_priority] > \
                        priority_order[tokens[2]]:
                    pass
                elif priority_order[highest_priority] == \
                        priority_order[tokens[2]]:
                    gene_list.append(tokens[3])
                else:  # higher priority
                    gene_list = [tokens[3]]
                    highest_priority = tokens[2]
        # for last peak
        writer.write("{}\t{}\t{}\n".format(
            current_peak,
            highest_priority,
            random.choice(gene_list)))
