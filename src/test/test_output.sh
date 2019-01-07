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

time python ../assign_peaks.py \
    -i "cat_R18_R25_all_groups_all_peaks_SORTED_MERGED_peak_centers.txt" \
    -a "test_center.txt" \
    -r "GCF_003254395.2_Amel_HAv3.1_genomic_PLUS_UTRs_introns.gff" \
    -p "test_priority.txt" \
    -z 1000

cmp test_center.txt peak_centers_all_new.txt || exit 1
echo center passed!

cmp test_priority.txt peak_centers_annotated_all_features_PRIORITY_ONLY.txt || exit 1

echo priority passed!
