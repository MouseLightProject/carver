#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:49:11 2019

@author: taylora
"""

import os
import numpy as np
import skimage.io as skio


def read_parameter_file(parameter_file_name):
    # reads calculated_parameters.txt file and parse it into a transform
    # const jobname = "ocJHfFH"
    # const nlevels = 6
    # const nchannels = 2
    # const shape_leaf_px = [406,256,152]
    # const voxelsize_used_um = [0.329714,0.342888,1.00128]
    # const origin_nm = [61677816,45421726,16585827]
    # const tile_type = convert(Cint,1)
    # const render_version = "2017-03-06 14:11:38 -0500 ec13bbfa7f9285447d3b9702b96a1f1afb847244"
    # const mltk_bary_version = "2016-11-11 12:16:03 -0500 84e153640047e3830abf835e1da4b738efa679d3"
    # const tilebase_version = "2016-08-22 15:49:39 -0400 cc171869a904e9e876426f2bb2732a38e607a102"
    # const nd_version = "2016-11-17 08:30:01 -0500 ef4923831c7bddadd0bba6b706f562a7cde00183"
    # const ndio_series_version = "2016-08-23 11:11:13 -0400 fdfe30a71f3d97fad6ac9982be50d8aea90b5234"
    # const ndio_tiff_version = "2016-08-23 11:11:54 -0400 df46d485cdf55ba66b8ed16fcf9fd9f3d5892464"
    # const ndio_hdf5_version = "2016-08-30 14:25:54 -0400 0c7ac77c5ca535913bfae5300159e6bdf60e36ca"
    # const mylib_version = "2013-08-06 19:15:35 -0400 0ca27aae55a5bab44263ad2e310e8f4507593ddc"
    params = {} # initialize dictionary
    with open(parameter_file_name, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            parts = text.split('=')
            keyval = parts[0].strip()
            if keyval == 'const nlevels':
                params['level_count'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const shape_leaf_px':
                params['leaf_shape'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const voxelsize_used_um':
                params['spacing'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const origin_nm':
                params['origin'] = 0.001 * np.array(eval(parts[1].strip('\n')),dtype=np.float)  # convert to um
            elif keyval == 'const nchannels':
                params['channel_count'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
    return params


def ijk2oct(ijk, level_count, chunk_shape):
    # converts ijk location to oct location
    # ALT: params['nlevels'] is the number of levels in the octree, i.e. the length of the path to a leaf stack
    # ALT: params['leafshape'] is the shape of the leaf stacks in the octree, in ijk/xyz order
    # ALT: ijk is in voxels, and is zero-based.  Values should be integral, although the type need not be, I don't think.
    # ALT" ijk is n x 3
    # ALT: chunk_shape is the chunk shape at the bottom level of the tree, i.e. the most zoomed-in level
    # ALT: Returns two things: First is the octree path, an n x nlevels array, row i the path to the leaf stack for ijk[i]
    # ALT:                     Second n x 3, row i giving the ijk position in the leaf stack for ijk[i]
    if ijk.ndim==1:
        did_promote = True
        ijk = ijk[None,:]  # convert vector to a single-row matrix
    else:
        did_promote = False

    level_count = np.int(level_count)
    #chunk_shape = params['leafshape']
    row_count = ijk.shape[0]
    octpath = np.zeros((row_count,level_count), dtype='int64')
    ijk_within_chunk = np.zeros((row_count,3), dtype='int64')
    for idx in range(row_count):
        bits = []
        ijk_this = ijk[idx]
        #u = chunk_shape
        for n in range(level_count-1,-1,-1):
            bn = (2**n)*chunk_shape
            th = ijk_this>bn
            bits.append(th)
            ijk_this = ijk_this - bn*th
        # convert to octodigit
        octpath[idx,:] = (1 + np.sum(np.array(bits)*2**np.array([0,1,2]),axis=1))[None,:]
        ijk_within_chunk[idx,:] = ijk_this

    if did_promote:
        # demote back to vectors
        octpath = np.ndarray.flatten(octpath)
        ijk_within_chunk = np.ndarray.flatten(ijk_within_chunk)

    return octpath, ijk_within_chunk


def load_mouselight_rendered_chunk(rendered_folder_name, chunk_path, channel_index):
    chunk_path_as_string = '/'.join(str(a) for a in chunk_path)
    chunk_folder_name = os.path.join(rendered_folder_name, chunk_path_as_string)
    chunk_file_leaf_name = 'default.%d.tif' % channel_index
    chunk_file_name = os.path.join(chunk_folder_name, chunk_file_leaf_name)
    raw_stack = skio.imread(chunk_file_name, plugin="tifffile")
    stack = np.transpose(raw_stack)
    return stack


def get_mouselight_rendered_substack(rendered_folder_name, channel_index, stack_origin, stack_shape):
    # stack_origin should be in zero-based voxel coords
    parameter_file_name = os.path.join(rendered_folder_name, 'calculated_parameters.jl')
    parameters = read_parameter_file(parameter_file_name)
    level_count = parameters['level_count']
    chunk_shape = parameters['leaf_shape']
    (chunk_path, stack_origin_chunk_ijk) = ijk2oct(stack_origin, level_count, chunk_shape)
    stack_far_corner_chunk_ijk = stack_origin_chunk_ijk + stack_shape
    if np.any(stack_far_corner_chunk_ijk-1 > chunk_shape):
        raise RuntimeError('Requested stack has to be contained within one chunk')
    chunk = load_mouselight_rendered_chunk(rendered_folder_name, chunk_path, channel_index)
    stack = chunk[stack_origin_chunk_ijk[0]:stack_far_corner_chunk_ijk[0],
                  stack_origin_chunk_ijk[1]:stack_far_corner_chunk_ijk[1],
                  stack_origin_chunk_ijk[2]:stack_far_corner_chunk_ijk[2]]
    return stack

                  

    

import matplotlib
import matplotlib.pyplot as plt
import z5py
import util

#matplotlib.use('Qt5Agg')

def real_zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing):
    return (sample_point_xyz - origin) / spacing


def zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing):
    ijk_real = real_zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing)
    return (np.round(ijk_real)).astype('int')  # must be round(), not floor()


octree_path = '/nrs/mouselight/SAMPLES/2018-10-01'

#n5_file_path = '/home/taylora/cache/gt/2018-10-01/navigator-output/test-swcs-carved.n5'
#n5_file_path = '/groups/scicompsoft/home/taylora/cache/gt/2018-10-01/navigator-output/test-swcs-carved.n5'
#offset = [75445.691517998261       16536.315191998805       36051.105973001562]   # soma 1
#offset = np.array([75885.516008994338,      15382.253859001321,       35947.93727700039])   # soma 2
#sample_point_xyz_centered = [2074.0422090017528      -684.41304799880709    1333.5478669984368] ;  # soma 1
#sample_point_xyz_centered = np.array([1157.5993180056539,       11.56324699867946,      1683.9443079996054])   # soma 2
#sample_point_xyz_centered = [738.19657500174071      -1424.7829449988058       702.41362699843739] ;

#n5_file_path = '/nrs/mouselight/cluster/navigator-output/2018-10-01-tile-chunks-on-cluster/consensus-neurons-with-machine-centerpoints-labelled-as-swcs-carved.n5'
#offset = np.array([75621.168080000323,       16211.625832000847,       36208.061435002361])   # G-040
#sample_point_xyz_centered = np.array([-501.96361500032071,      -1013.1980690008477,      -4455.3761980023555])   # G-040 centerpoint 33650

#n5_file_path = '/nrs/mouselight/cluster/navigator-output/2018-10-01-tile-chunks-on-cluster/consensus-neurons-with-machine-centerpoints-labelled-as-swcs-carved.n5'
#offset = np.array([0, 0 ,0])
#sample_point_xyz_centered = np.array([75148.9, 15244.3, 31761.2])   # G-040 centerpoint copied from WS

n5_file_path = '/nrs/mouselight/cluster/navigator-output/2018-10-01-tile-chunks-on-cluster/consensus-neurons-with-machine-centerpoints-labelled-as-swcs-carved.n5'
#offset = np.array([73814.192263002085,       17366.654304000425,       36205.867519004452])
#sample_point_xyz_centered = np.array([479.90134299792408,        436.8809699995727,      -1025.1585660044511])   # G-017 centerpoint 6320

voxel_center_in_swc = np.array([74447.743653, 17765.863807, 34985.740732])
# this point survives import and export unchanged, and is shown as a voxel center in x-y
# it is also s.t. (it - origin_jaws)/spacing is very nearly integral


#[origin, spacing] = load_transform_txt(fullfile(octree_path, 'transform.txt')) ;
params = util.readParameterFile(parameterfile = os.path.join(octree_path, 'calculated_parameters.jl'))
tile_level_count = params["nlevels"].astype(int)
tile_shape = params["leafSize"].astype(int)
#origin = params["origin"]
spacing = params["spacing"]
origin_jaws = np.array([69445.01944245792, 12917.29486937314, 30198.96941474185])
origin = origin_jaws

sample_point_xyz = voxel_center_in_swc  # this is a voxel center


channel_index = 0 
infinity_norm_radius = 4  # we want a cube of diameter twice this, in radial voxel units
real_sample_point_ijk = real_zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing)
sample_point_ijk = zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing)
#half_diagonal = [512 512 128] ;
half_diagonal = (infinity_norm_radius * np.array([1, 1, 1/4])).astype('int')
stack_lower_corner = sample_point_ijk - half_diagonal
stack_shape = 2*half_diagonal
stack_upper_corner = stack_lower_corner + stack_shape
  # Since all elements of stack_shape are even, our pixel cannot strictly be at the center of the image.
  # This puts our voxel center one half-pixel to the left, and one half-pixel down, and one half-pixel 'away' from the
  # center of the stack.  In particular, it should be at stack[infinity_norm_radius, infinity_norm_radius, infinity_norm_radius/4].

f = z5py.File(n5_file_path, 'r')
dset = f['volume']
shape = dset.shape
dtype = dset.dtype

raw_stack_4d = dset[stack_lower_corner[0]:stack_upper_corner[0],
                    stack_lower_corner[1]:stack_upper_corner[1],
                    stack_lower_corner[2]:stack_upper_corner[2],
                    channel_index]
stack = np.squeeze(raw_stack_4d, axis=3)

plt.close('all')


slice = stack[:,:,int(infinity_norm_radius/4)]
plt.figure()
plt.imshow(slice.T)
plt.show()

# Now get data direct from the slice
direct_stack = get_mouselight_rendered_substack(octree_path, channel_index, stack_lower_corner, stack_shape)
direct_slice = direct_stack[:,:,int(infinity_norm_radius/4)]
plt.figure()
plt.imshow(direct_slice.T)

plt.show()

assert(np.all(np.ndarray.flatten(np.equal(slice, direct_slice))))



    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        