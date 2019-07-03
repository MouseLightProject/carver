import improc
import os
import util
import numpy as np
import pickle

def crop_from_render(render_folder_name, input_swc_file_or_folder_name, output_folder_name, output_volume_file_name, do_use_simple_for_loop=False):
    output_volume_file_path =  os.path.join(output_folder_name, output_volume_file_name)

    params = util.readParameterFile(parameterfile=render_folder_name + "/calculated_parameters.jl")
    tile_level_count = params["nlevels"].astype(int)
    tile_shape = params["leafSize"].astype(int)
    origin_um = params["origin"]
    spacing_um = params["spacing"]
    full_volume_shape = tile_shape * 2**tile_level_count

    # check if input argument is file or folder
    if os.path.isfile(input_swc_file_or_folder_name):
        inputfolder, swc_name_w_ext = os.path.split(input_swc_file_or_folder_name)
        xyz_um, edges, R, offset, scale, header = util.readSWC(os.path.join(inputfolder, swc_name_w_ext))
    elif os.path.isdir(input_swc_file_or_folder_name):
        inputfolder = input_swc_file_or_folder_name
        xyz_um, edges, R = util.appendSWCfolder(inputfolder) # somewhat redundant but cleaner
        xyz_um_, edges_, R_, filenames, header = util.readSWCfolder(inputfolder)
    else:
        raise RuntimeError('%s does not seem to be a file, nor a folder' % input_swc_file_or_folder_name)

    # Convert swc coords to voxels
    xyz_in_voxels = util.um2pix(xyz_um, origin_um, spacing_um)

    # Each tile in the rendered image will itself be 'octreed' into a set of 'leafs'
    # depthextend tells now many octree levels there will be within each tile
    extra_level_count = 3
    leaf_level_count = tile_level_count + extra_level_count
    leaf_shape = (tile_shape / (2**extra_level_count)).astype(int)
    octpath, xres = improc.ijk2oct(xyz_in_voxels, leaf_level_count, leaf_shape)

    #depthFull = params_p1["nlevels"].astype(int)
    #leaf_shape = params_p1["leafshape"].astype(int)

    swc_base_name = os.path.basename(input_swc_file_or_folder_name)
    tile_list_pickle_file_name = '%s-tile-list.pickle' % swc_base_name
    tile_list_pickle_file_path = os.path.join(output_folder_name, tile_list_pickle_file_name)
    try:
        tile_hash = pickle.load(open(tile_list_pickle_file_path, 'rb'))
        print('Loaded tile list from memo file')
    except os.error:
        octpath_cover = np.unique(octpath, axis=0)
        gridlist_cover = improc.oct2grid(octpath_cover)

        print('About to start dilation...')
        #octpath_dilated = octpath_cover.copy()
        desired_carve_out_half_diagonal_as_scalar = 512
        desired_carve_out_half_diagonal = desired_carve_out_half_diagonal_as_scalar * np.array([1.0, 1.0, 1.0/4.0])
        #dilation_count = 8
        dilation_count = np.max( np.ceil(desired_carve_out_half_diagonal.astype(float) / leaf_shape.astype(float)) ).astype(int).item()
        # should be enough to get about a 512 vx cube around each swc centerpoint
        # (except 4x less in z, b/c axial rez is less)
        octpath_dilated, junk = improc.dilateOct(octpath_cover, dilation_count)
        #for dilation_index in range(dilation_count):
        #    print('Finished dilation iteration %d of %d' % (dilation_index+1, dilation_count))
        print('Done with dilation!')

        tile_hash = improc.chunklist(octpath_dilated, tile_level_count) #1..8
        pickle.dump(tile_hash, open(tile_list_pickle_file_path, 'wb'))

    #tileids = list(tile_hash.keys())
    # base on bounding box (results in cropped output volume)
    # gridReference = np.min(gridlist_dilated, axis=0)
    # gridSize = np.max(gridlist_dilated, axis=0) - gridReference +1
    # base on initial image
    #gridReference = np.array((0,0,0))
    #gridSize = full_volume_shape/leaf_shape
    #   # 3-array, number of leaves in each dimension to make up the full volume

    #volReference = gridReference*leaf_shape
    #full_volume_shape_including_color_channel = np.append(full_volume_shape, 2)  # append color channel
    #chunk_shape_including_color_channel = np.append(leaf_shape, 2)
    #chunk_shape_including_color_channel = np.append(tile_shape, 2)

    # setting = dict()
    # setting['volSize'] = full_volume_shape
    # setting['chunkSize'] = chunksize
    # setting['depthBase'] = tile_level_count
    # setting['depthFull'] = leaf_level_count
    # setting['tileSize'] = tile_shape
    # setting['leaf_shape'] = leaf_shape
    #
    # setting['dtype'] = 'uint16'

    output_file_extension = os.path.splitext(output_volume_file_name)[1]
    if output_file_extension == '.h5' :
        output_file_type = 'h5'
        compression_method = "gzip"
        compression_options = 9
    elif output_file_extension == '.n5' :
        output_file_type = 'n5'
        compression_method = "gzip"
        compression_options = {'level': 9}
    elif output_file_extension == '.zarr':
        output_file_type = 'zarr'
        compression_method = "blosc"
        compression_options = {}
    else :
        raise RuntimeError('Don''t recognize the output file extension %s' % output_file_extension)

    # Finally, write the voxel carved data to disk
    color_channel_count = 2
    util.dump_write(render_folder_name,
                    full_volume_shape,
                    'uint16',
                    color_channel_count,
                    output_volume_file_path,
                    tile_hash,
                    leaf_level_count,
                    tile_level_count,
                    compression_method,
                    compression_options,
                    output_file_type,
                    do_use_simple_for_loop)

# end def crop_from_render()
