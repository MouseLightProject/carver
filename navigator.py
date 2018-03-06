# -*- coding: utf-8 -*-
"""Command line interface for navigator."""

import util
import numpy as np
import cropper
import os
import sys, getopt
import stat
import h5py

import shutil

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def fixKinksinAnnotation():
    input_folder = '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo'
    swcfiles = [os.path.join(input_folder, fold, files) for fold in os.listdir(input_folder) if
                os.path.isdir(os.path.join(input_folder, fold)) for files in
                os.listdir(os.path.join(input_folder, fold)) if
                files.endswith("-carved.swc")]
    swcfiles.sort()
    for swc_file in swcfiles[1]:
        path, filename = os.path.split(swc_file)
        output_h5_file = os.path.join(path, filename.split('.')[0] + '.h5')
        input_swc = swc_file
        # output_h5_file = os.path.join(path,filename.split('.')[0][:-1]+'.h5')
        # input_swc = os.path.join(path,filename.split('.')[0][:-1]+'.swc')
        with h5py.File(output_h5_file, "r+") as f:
            try:
                del f['reconstruction']
            except Exception:
                pass
            um, edges, R, offset, scale, header = util.readSWC(swcfile=input_swc, scale=1)
            dset_swc = f.create_dataset("reconstruction", (um.shape[0], 7), dtype='f')
            for iter, xyz_ in enumerate(um):
                xyz_ = np.ceil(xyz_ - np.sqrt(np.finfo(float).eps))
                dset_swc[iter, :] = np.array(
                    [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])


def main(argv):
    # 'python navigator.py -i/nrs/mouselight/SAMPLES/2017-11-17 -s /groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/recon_repo/2017-11-17_G-017_Seg-3.swc -o /groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data/2017-11-17_G-017_Seg-3'
    # JW octree depth
    number_of_level = 3
    try:
        opts, args = getopt.getopt(argv,"hi:s:o:",["data_fold=","input_swc_file=","output_folder="])
    except getopt.GetoptError:
        print('navigator.py -i <data_folder> -s <swc_file> -o <output_folder>')
        sys.exit(2)
    data_fold='/nrs/mouselight/SAMPLES/2017-09-25-padded'
    input_swc_file='/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo/2017-09-25_G-001_consensus/2017-09-25_G-001_consensus-carved.swc'
    output_folder='/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo/2017-09-25_G-001_consensus'
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

    inputfolder, swc_name_w_ext = os.path.split(input_swc_file)
    swc_name, _ = swc_name_w_ext.split(os.extsep)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    output_swc_name = '{}-carved.swc'.format(swc_name)
    output_h5_name =  '{}-carved.h5'.format(swc_name)
    JW_output_folder = os.path.join(output_folder,'JW')
    if not os.path.exists(JW_output_folder):
        cropper.crop_from_render(data_fold, input_swc_file, output_folder, output_swc_name, output_h5_name)
    # shutil.rmtree(JW_output_folder)
    if not os.path.exists(JW_output_folder):
        os.makedirs(JW_output_folder)
        os.chmod(JW_output_folder, 0o770)

    output_h5_file = os.path.join(output_folder, output_h5_name)
    converter = util.Convert2JW(output_h5_file, JW_output_folder, number_of_level=None)
    converter.convert2JW()
    converter.mergeJW(number_of_level=converter.number_of_level)
    converter.create_transform_file()
    print('DONE')

if __name__ == "__main__":
   main(sys.argv[1:])

# [21118.5, 4641.6, 3385.0]
#
#
#
#
#
