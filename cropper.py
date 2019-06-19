import improc
import os
import util
import re
import numpy as np
import h5py
import skimage.io as io

def crop_from_render(data_fold, input_swc,output_folder, output_swc_name, output_h5_name):
    output_swc_file = os.path.join(output_folder,output_swc_name)
    output_h5_file =  os.path.join(output_folder,output_h5_name)

    params = util.readParameterFile(parameterfile=data_fold+"/calculated_parameters.jl")
    tile_level_count = params["nlevels"].astype(int)
    tile_shape = params["leafSize"].astype(int)
    origin_um = params["origin"]
    spacing_um = params["spacing"]

    # check if input argument is file or folder
    if os.path.isfile(input_swc):
        inputfolder, swc_name_w_ext = os.path.split(input_swc)
        xyz_um, edges, R, offset, scale, header = util.readSWC(os.path.join(inputfolder, swc_name_w_ext))
    elif os.path.isdir(input_swc):
        inputfolder = input_swc
        xyz_um, edges, R = util.appendSWCfolder(inputfolder) # somewhat redundant but cleaner
        xyz_um_, edges_, R_, filenames, header = util.readSWCfolder(inputfolder)
    else:
        raise RuntimeError('%s does not seem to be a file, nor a folder' % input_swc)

    # Convert swc coords to voxels
    xyz_in_voxels = util.um2pix(xyz_um, origin_um, spacing_um)

    # Each tile in the rendered image will itself be 'octreed' into a set of 'leafs'
    # depthextend tells now many octree levels there will be within each tile
    extra_level_count = 3
    leaf_level_count = tile_level_count + extra_level_count
    leaf_shape = (tile_shape / (2**extra_level_count)).astype(int)
    octpath, xres = improc.xyz2oct(xyz_in_voxels, leaf_level_count, leaf_shape)

    #depthFull = params_p1["nlevels"].astype(int)
    #leaf_shape = params_p1["leafshape"].astype(int)

    octpath_cover = np.unique(octpath, axis=0)
    gridlist_cover = improc.oct2grid(octpath_cover)

    octpath_dilated = octpath_cover.copy()
    desired_carve_out_half_diagonal_as_scalar = 512
    desired_carve_out_half_diagonal = desired_carve_out_half_diagonal_as_scalar * np.array([1.0, 1.0, 1.0/4.0])
    #dilation_count = 8
    dilation_count = np.max( np.ceil(desired_carve_out_half_diagonal.astype(float) / leaf_shape.astype(float)) ).astype(int).item()
    # should be enough to get about a 512 vx cube around each swc centerpoint
    # (except 4x less in z, b/c axial rez is less)
    for dilation_index in range(dilation_count):
        octpath_dilated, gridlist_dilated = improc.dilateOct(octpath_dilated)


    tilelist = improc.chunklist(octpath_dilated, tile_level_count) #1..8

    tileids = list(tilelist.keys())
    # base on bounding box (results in cropped output volume)
    # gridReference = np.min(gridlist_dilated, axis=0)
    # gridSize = np.max(gridlist_dilated, axis=0) - gridReference +1
    # base on initial image
    gridReference = np.array((0,0,0))
    gridSize = tile_shape*(2**(tile_level_count))/leaf_shape

    volReference = gridReference*leaf_shape
    outVolumeSize = np.append(gridSize*leaf_shape,2) #append color channel
    chunksize = np.append(leaf_shape,2)

    setting = dict()
    setting['volSize'] = outVolumeSize
    setting['chunkSize'] = tuple(chunksize)
    setting['depthBase'] = tile_level_count
    setting['depthFull'] = leaf_level_count
    setting['tileSize'] = tile_shape
    setting['leaf_shape'] = leaf_shape
    setting['volReference'] = volReference

    setting['compression'] = "gzip"
    setting['compression_opts'] = 9
    setting['dtype'] = 'uint16'
    setting['type'] = 'h5'

    # Finally, write the voxel carved data to disk
    util.dump_write(data_fold, output_h5_file, setting, tilelist)

# end def crop_from_render()
