import util
import numpy as np
import matplotlib.pyplot as plt
import improc
import cropper
import os

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

if __name__ == '__main__':

    # oct in [1...8]
    # grid in [0...(2**depth-1)]

    import importlib
    importlib.reload(util)
    # importlib.reload(improc)

    data_fold = '/nrs/mouselight/SAMPLES/2017-11-17'
    swc_folder = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/recon_repo'
    swc_name = '2017-11-17_G-017_Seg-1'
    input_swc_file =  os.path.join(swc_folder,swc_name+'.swc')

    output_folder = os.path.join('/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data',swc_name)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    output_swc_name = '{}-cropped.swc'.format(swc_name)
    output_h5_name =  '{}-cropped.h5'.format(swc_name)
    cropper.crop_from_render(data_fold,input_swc_file,output_folder,output_swc_name,output_h5_name)

    JW_output_folder = os.path.join(output_folder,'JW')
    if not os.path.exists(JW_output_folder):
        os.makedirs(JW_output_folder)
    number_of_level = 3
    output_h5_file = os.path.join(output_folder,output_h5_name)
    converter = util.Convert2JW(output_h5_file,JW_output_folder,number_of_level)
    converter.convert2JW()
    converter.mergeJW(number_of_level=3)
    converter.create_transform_file()
