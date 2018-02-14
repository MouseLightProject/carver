# -*- coding: utf-8 -*-
"""Command line interface for navigator."""

import util
import numpy as np
import cropper
import os
import sys, getopt

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def main(argv):
    # 'python navigator.py -i/nrs/mouselight/SAMPLES/2017-11-17 -s /groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/recon_repo/2017-11-17_G-017_Seg-3.swc -o /groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data/2017-11-17_G-017_Seg-3'
    # JW octree depth
    number_of_level = 3
    try:
        opts, args = getopt.getopt(argv,"hi:s:o:",["data_fold=","input_swc_file=","output_folder="])
    except getopt.GetoptError:
        print('navigator.py -i <data_folder> -s <swc_file> -o <output_folder>')
        sys.exit(2)
    input_swc_file=''
    output_folder=''
    data_fold=''
    for opt, arg in opts:
        print('opt:', opt,'arg:', arg)
        if opt == '-h':
            print('navigator.py -i <data_folder> -s <swc_file> -o <output_folder>')
            sys.exit()
        elif opt in ("-i", "--data_fold"):
            print(arg)
            data_fold = arg
        elif opt in ("-s", "--input_swc_file"):
            input_swc_file = arg
        elif opt in ("-o", "--output_folder"):
            output_folder = arg

    print('SWCFILE   :', input_swc_file)
    print('DATAFOLDER   :', data_fold)
    print('OUTPUT    :', output_folder)

    # oct in [1...8]
    # grid in [0...(2**depth-1)]

    inputfolder,swc_name_w_ext = os.path.split(input_swc_file)
    swc_name,_ = swc_name_w_ext.split(os.extsep)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    output_swc_name = '{}-carved.swc'.format(swc_name)
    output_h5_name =  '{}-carved.h5'.format(swc_name)
    cropper.crop_from_render(data_fold,input_swc_file,output_folder,output_swc_name,output_h5_name)

    JW_output_folder = os.path.join(output_folder,'JW')
    if not os.path.exists(JW_output_folder):
        os.makedirs(JW_output_folder)
    output_h5_file = os.path.join(output_folder,output_h5_name)
    converter = util.Convert2JW(output_h5_file,JW_output_folder,number_of_level)
    converter.convert2JW()
    converter.mergeJW(number_of_level=number_of_level)
    converter.create_transform_file()

if __name__ == "__main__":
   main(sys.argv[1:])
