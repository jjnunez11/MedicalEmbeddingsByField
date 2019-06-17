from __future__ import division
import argparse
from analysis_choi_mrp import print_choi_mrp
from analysis_choi_mcsp import print_choi_mcsp
from analysis_yu_umnsrs_cor import print_yu_umnsrs_cor
from analysis_beam_bootstrap import print_beam_bootstrap
from analysis_new_sysvec import print_new_sysvec

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filenames", default='config_embed_files.txt')
    args = parser.parse_args()
    filenames_file = args.filenames

    filenames = []
    with open(filenames_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            filenames.append(row.strip().split(','))
    
    print_choi_mrp(filenames) 
    #print_yu_umnsrs_cor(filenames)
    #print_choi_mcsp(filenames) 
    #print_beam_bootstrap(filenames) 
    #print_new_sysvec(filenames)
