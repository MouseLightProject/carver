import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import z5py
import util

matplotlib.use('Qt5Agg')

def zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing):
    return (np.floor((sample_point_xyz - origin) / spacing)).astype('int')

octree_path = '/nrs/mouselight/SAMPLES/2018-07-02'

#n5_file_path = '/home/taylora/cache/gt/2018-10-01/navigator-output/test-swcs-carved.n5'
#n5_file_path = '/groups/scicompsoft/home/taylora/cache/gt/2018-10-01/navigator-output/test-swcs-carved.n5'
#offset = [75445.691517998261       16536.315191998805       36051.105973001562]   # soma 1
#offset = np.array([75885.516008994338,      15382.253859001321,       35947.93727700039])   # soma 2
#sample_point_xyz_centered = [2074.0422090017528      -684.41304799880709    1333.5478669984368] ;  # soma 1
#sample_point_xyz_centered = np.array([1157.5993180056539,       11.56324699867946,      1683.9443079996054])   # soma 2
#sample_point_xyz_centered = [738.19657500174071      -1424.7829449988058       702.41362699843739] ;

n5_file_path = '/groups/mousebrainmicro/mousebrainmicro/scripts/carver/2018-07-02-test-swcs-output/fluorescence-near-consensus-parallel.n5'
#offset = np.array([76542.556325, 18907.477255, 33767.499672])   # G-007
#sample_point_xyz_centered = np.array([76174.5, 18786.1, 34071.6])   # G-007 node near an end point

# offset the sample point
sample_point_xyz = np.array([76174.5, 18786.1, 34071.6])   # G-007 node near an end point

#[origin, spacing] = load_transform_txt(fullfile(octree_path, 'transform.txt')) ;
params = util.readParameterFile(parameterfile = os.path.join(octree_path, 'calculated_parameters.jl'))
tile_level_count = params["nlevels"].astype(int)
tile_shape = params["leafSize"].astype(int)
origin = params["origin"]
spacing = params["spacing"]

sample_point_ijk = zero_based_ijk_from_xyz(sample_point_xyz, origin, spacing)
#half_diagonal = [512 512 128] ;
half_diagonal = (1024 * np.array([1, 1, 1/4])).astype('int')
stack_lower_corner = sample_point_ijk - half_diagonal
stack_shape = 2*half_diagonal
stack_upper_corner = stack_lower_corner + stack_shape

f = z5py.File(n5_file_path, 'r')
dset = f['volume']
shape = dset.shape  # this will be czyx order, b/c n5
dtype = dset.dtype

raw_stack_4d = dset[0,
                    stack_lower_corner[2]:stack_upper_corner[2],
                    stack_lower_corner[1]:stack_upper_corner[1],
                    stack_lower_corner[0]:stack_upper_corner[0]]
stack_transposed = np.squeeze(raw_stack_4d, axis=0)
stack = np.transpose(stack_transposed)

slice = stack[:,:,256]
imgplot = plt.imshow(slice.T)

(m,n) = slice.shape
central_radius = 100
central_slice = slice[int(m/2)-central_radius:int(m/2)+central_radius, int(n/2)-central_radius:int(n/2)+central_radius]
imgplot2 = plt.imshow(central_slice.T)

plt.show()
