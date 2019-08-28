#!/bin/bash

#self.infile = "cat_R18_R25_all_groups_all_peaks_SORTED_MERGED_peak_centers.txt"
#self.outfile = "peak_centers_annotated_all_features.txt"
#self.ref_gff = "GCF_003254395.2_Amel_HAv3.1_genomic_PLUS_UTRs_introns.gff"
#self.promoter_size = 1000
#
#infile = r'peak_centers_annotated_all_features.txt'
#outfile = r'peak_centers_annotated_all_features_PRIORITY_ONLY_try2.txt'

set -euo pipefail

[[ -s test_center.txt ]] && rm test_center.txt
[[ -s test_priority.txt ]] && rm test_priority.txt

python ../sorted_assign_gtf.py \
    -g dmel-all-r6.29-sort.gtf \
    -l arnold_S2_enhancers_peak_centers_v6_sort.txt \
    -o test_center.txt \
    -p test_priority.txt \
