#!/bin/bash

#self.infile = "cat_R18_R25_all_groups_all_peaks_SORTED_MERGED_peak_centers.txt"
#self.outfile = "peak_centers_annotated_all_features.txt"
#self.ref_gff = "GCF_003254395.2_Amel_HAv3.1_genomic_PLUS_UTRs_introns.gff"
#self.promoter_size = 1000
#
#infile = r'peak_centers_annotated_all_features.txt'
#outfile = r'peak_centers_annotated_all_features_PRIORITY_ONLY_try2.txt'

time python ../USE_FIRST_assign_peaks_all_features.py \
    -i "cat_R18_R25_all_groups_all_peaks_SORTED_MERGED_peak_centers.txt" \
    -o "test_center.txt" \
    -r "GCF_003254395.2_Amel_HAv3.1_genomic_PLUS_UTRs_introns.gff" \
    -p 1000

time python ../USE_SECOND_pull_out_priority_features.py \
    -i test_center.txt \
    -o test_priority.txt

cmp test_priority.txt peak_centers_annotated_all_features_PRIORITY_ONLY.txt || exit 1

echo passed!
