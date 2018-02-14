# -*- coding: utf-8 -*-
"""Creates a shell script to run batch process on command line."""

import os

def main(argv):
    # parses input repo to generate list of consensus swcs
    function_path = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/navigator.py'
    input_folder = '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/swcrepo/2017-09-25'
    sample_folder = '/nrs/mouselight/SAMPLES/2017-09-25-padded'
    output_folder = '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo'
    output_sh_file = '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo/run_navigator.sh'

    consensus_swc_files = [os.path.join(root, name) for root, dirs, files in os.walk(input_folder) for name in files if name.endswith(("Consensus.swc", "consensus.swc"))]
    consensus_swc_files.sort()

    # for each consensus file run navigator script
    # usage:
    # print('navigator.py -i <data_folder> -s <swc_file> -o <output_folder>')

    with open(output_sh_file,'w') as fswc:
        for swc_file in consensus_swc_files:
            path, filename = os.path.split(swc_file)
            mystr = 'python {} -i {} -s {} -o {}&\n'.format(function_path,sample_folder, swc_file,
                                                              os.path.join(output_folder, filename.split('.')[0]))
            fswc.write(mystr)

    os.system('chmod g+x '+output_sh_file)
