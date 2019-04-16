from __future__ import division
import argparse
from analysis_choi_mrp import print_choi_mrp
from analysis_choi_mcsp import print_choi_mcsp
from analysis_yu_umnsrs_cor import print_yu_umnsrs_cor
from analysis_beam_bootstrap import print_beam_bootstrap
from analysis_new_centroid import print_new_centroid

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filenames", default='orig_files_all.txt')
    args = parser.parse_args()
    filenames_file = args.filenames

    filenames = []
    with open(filenames_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            filenames.append(row.strip().split(','))
    
    
    
    ##filename_to_print, finegrain_ndcgs = get_css_analysis(filenames, num_of_nn, 'f')
    ##filename_to_print, coarsegrain_ndcgs = get_css_analysis(filenames, num_of_nn, 'c')
    ##for name, finegrain, coarsegrain in zip(filename_to_print, finegrain_ndcgs, coarsegrain_ndcgs):
    ##    print '%s & %.2f & %.2f \\\\' %(name.split('/')[-1], finegrain*100, coarsegrain*100)
    
    #print_choi_mrp(filenames) CHECKED
    print_yu_umnsrs_cor(filenames)
    #print_choi_mcsp(filenames) CHECKED
    #print_beam_bootstrap(filenames) CHECKED
    #print_new_centroid(filenames)
    
    
    ##filename_to_print, psych_ndcgs = get_choi_mrp_by_system(filenames, num_of_nn, 'c')
    ##filename_to_print, coarsegrain_ndcgs = get_css_analysis(filenames, num_of_nn, 'c')
    ##for name, psychs, coarsegrain in zip(filename_to_print, psych_ndcgs, coarsegrain_ndcgs):
    ##    print '%s & %.2f & %.2f \\\\' %(name.split('/')[-1], psychs*100, coarsegrain*100)
