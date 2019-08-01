# -*- coding: utf-8 -*-
"""Command line interface for navigator."""

import util
import numpy as np
import os
import sys
import getopt
import h5py
import improc
import pickle
from functools import partial
import z5py
import tqdm
from dask_jobqueue import LSFCluster
from dask.distributed import Client, LocalCluster
from dask.distributed import wait
from dask.distributed import progress
import getpass
import time



def dump_single_tile_id(tile_id, 
                        leaf_ids_within_tile, 
                        rendered_folder_path, 
                        tile_shape, 
                        leaf_shape, 
                        chunk_shape_with_color_as_tuple, 
                        dtype, 
                        dataset, 
                        is_dataset_transposed):
    tilename = '/'.join(a for a in tile_id)
    tilepath = os.path.join(rendered_folder_path, tilename)

    tile_octree_path = np.array(list(tile_id), dtype=int)
    tile_ijk_in_tile_grid = np.ndarray.flatten(improc.oct2grid(tile_octree_path.reshape(1, len(tile_octree_path))))
    tile_origin_ijk = tile_ijk_in_tile_grid * tile_shape
    tile_end_ijk = tile_origin_ijk + tile_shape
    if os.path.isdir(tilepath):
        tile_ijk_in_tile_grid_as_tuple = tuple(tile_ijk_in_tile_grid)
        tile_ijk_in_tile_grid_as_tuple_with_color = tile_ijk_in_tile_grid_as_tuple + (0,)  # dataset includes color channels
        if is_dataset_transposed:
            tile_ijk_in_tile_grid_as_tuple_with_color_maybe_flipped = tuple(reversed(tile_ijk_in_tile_grid_as_tuple_with_color)) 
        else:
            tile_ijk_in_tile_grid_as_tuple_with_color_maybe_flipped = tile_ijk_in_tile_grid_as_tuple_with_color
        does_chunk_exist = dataset.chunk_exists(tile_ijk_in_tile_grid_as_tuple_with_color_maybe_flipped)
        if not does_chunk_exist:
            im = improc.loadTiles(tilepath)
            # relativeDepth = leaf_level_count - tile_level_count
            output_tile_stack = np.zeros(chunk_shape_with_color_as_tuple, dtype=dtype)

            # patches in idTiled
            for leaf_octree_path_within_tile_as_string in leaf_ids_within_tile:
                leaf_octree_path_within_tile = np.array(list(leaf_octree_path_within_tile_as_string), dtype=int)
                leaf_ijk_in_leaf_grid_within_tile = improc.oct2grid(leaf_octree_path_within_tile.reshape(1, len(leaf_octree_path_within_tile)))  # in 0 base

                start = np.ndarray.flatten(leaf_ijk_in_leaf_grid_within_tile * leaf_shape)
                end = np.ndarray.flatten(start + leaf_shape)
                leaf_stack = im[start[0]:end[0], start[1]:end[1], start[2]:end[2], :]
                output_tile_stack[start[0]:end[0], start[1]:end[1], start[2]:end[2], :] = leaf_stack

            if is_dataset_transposed :
                dataset[:, tile_origin_ijk[2]:tile_end_ijk[2], tile_origin_ijk[1]:tile_end_ijk[1], tile_origin_ijk[0]:tile_end_ijk[0]] = \
                    np.transpose(output_tile_stack)                
            else:
                dataset[tile_origin_ijk[0]:tile_end_ijk[0], tile_origin_ijk[1]:tile_end_ijk[1], tile_origin_ijk[2]:tile_end_ijk[2], :] = \
                    output_tile_stack
# end def



def dump_write(render_folder_name,
               full_volume_shape,
               dtype,
               color_channel_count,
               output_file_name,
               tile_hash,
               leaf_level_count,
               tile_level_count,
               compression_method,
               compression_options,
               output_file_type,
               do_use_simple_for_loop=False):
    # dumps volumetric data into h5/n5/zarr
    #self.inputLoc = inputloc

    tile_shape = (full_volume_shape / (2**tile_level_count)).astype(int)
    leaf_shape = (full_volume_shape / (2**leaf_level_count)).astype(int)

    # check if dataset name is provided
    splitted_name = output_file_name.split(':')
    if  len(splitted_name) == 1:
        output_file_name =  splitted_name[0]
        dataset_name =  "volume"
    elif len(splitted_name) ==2:
        output_file_name =  splitted_name[0]
        dataset_name =  splitted_name[1]
    else:
        raise ValueError('output file name has more than one ":"', output_file_name)
    #self.setting = setting
    #self.tilelist = tilelist
    tile_id_list = list(tile_hash.keys())
    leaf_ids_per_tile_list = list(tile_hash.values())

    # # Unpack the settings
    # volSize = tuple(map(int,setting['volSize']))
    # tileSize = setting['tileSize']
    # #volReference = setting['volReference']
    # depthFull = setting['depthFull']
    # depthBase = setting['depthBase']
    # leafSize = setting['leaf_shape']
    # dtype = setting['dtype']
    # chunkSize = tuple(map(int,setting['chunkSize']))
    # compression_method = setting['compression']
    # comp_opts = setting['compression_opts']
    chunk_shape = tile_shape
    full_volume_shape_including_color_channel = np.append(full_volume_shape, color_channel_count)  # append color channel
    chunk_shape_including_color_channel = np.append(chunk_shape, color_channel_count)
    full_volume_shape_with_color_channels_as_tuple = tuple(map(int, full_volume_shape_including_color_channel))
    chunk_shape_with_color_as_tuple = tuple(map(int, chunk_shape_including_color_channel))

    if output_file_type=='h5':
        # write into h5
        with h5py.File(output_file_name, "w") as f:
            # dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0], 7), dtype='f')
            # for iter, xyz_ in enumerate(xyz_shifted):
            #     dset_swc[iter, :] = np.array(
            #         [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])
            dataset = f.create_dataset(dataset_name,
                                       full_volume_shape_with_color_channels_as_tuple,
                                       dtype=dtype,
                                       chunks=chunk_shape_with_color_as_tuple,
                                       compression=compression_method,
                                       compression_opts=compression_options)


            # crop chuncks from a tile read in tilelist
            for iter, tile_id in enumerate(tile_id_list):
                print('{} : {} out of {}'.format(tile_id, iter+1, len(tile_id_list)))
                leaf_id_within_tile = tile_hash[tile_id]
                dump_single_tile_id(tile_id,
                                    leaf_id_within_tile,
                                    render_folder_name,
                                    tile_shape,
                                    leaf_shape,
                                    chunk_shape_with_color_as_tuple,
                                    dtype,
                                    dataset, 
                                    is_dataset_transposed=False)
    elif output_file_type=='n5' or output_file_type=='zarr':
        # write into z5 or n5
        if do_use_simple_for_loop:
            use_zarr_format = (output_file_type == 'zarr')
            with z5py.File(output_file_name, 'a', use_zarr_format=use_zarr_format) as f:
                # require_dataset seems to choke on the compression_options {level: 9}, so this is a workaround
                g = f.require_group('/')
                try:
                    dataset = g[dataset_name]
                except KeyError:
                    dataset = f.create_dataset(dataset_name,
                                               shape=tuple(reversed(full_volume_shape_with_color_channels_as_tuple)),
                                               dtype=dtype,
                                               chunks=tuple(reversed(chunk_shape_with_color_as_tuple)),
                                               compression=compression_method,
                                               **compression_options)                    
                for tile_id in tqdm.tqdm(tile_id_list):
                    leaf_ids_within_tile = tile_hash[tile_id]
                    dump_single_tile_id(tile_id,
                                        leaf_ids_within_tile,
                                        render_folder_name,
                                        tile_shape,
                                        leaf_shape,
                                        chunk_shape_with_color_as_tuple,
                                        dtype,
                                        dataset,
                                        is_dataset_transposed=True)
        else:
            username = getpass.getuser()
            scratch_folder_path = '/scratch/%s' % username
            with LSFCluster(cores=1, memory='15 GB', local_dir=scratch_folder_path, projectstr='mouselight', queue='normal', extralist='-o /dev/null -e /dev/null') as cluster:
                cluster.adapt(minimum=1, maximum=1000)
                #cluster = LocalCluster(n_workers=4, threads_per_worker=1)
                #cluster.scale(200)
                with Client(cluster) as client:
                    use_zarr_format = (output_file_type=='zarr')
                    with z5py.File(output_file_name, 'a', use_zarr_format=use_zarr_format) as f:
                        # require_dataset seems to choke on the compression_options {level: 9}, so this is a workaround
                        g = f.require_group('/')
                        try:
                            dataset = g[dataset_name]
                        except KeyError:                        
                            dataset = f.create_dataset(dataset_name,
                                                       shape=tuple(reversed(full_volume_shape_with_color_channels_as_tuple)),
                                                       dtype=dtype,
                                                       chunks=tuple(reversed(chunk_shape_with_color_as_tuple)),
                                                       compression=compression_method,
                                                       **compression_options)
                        two_arg_dump_single_tile_id = \
                            partial(dump_single_tile_id,
                                    rendered_folder_path=render_folder_name,
                                    tile_shape=tile_shape,
                                    leaf_shape=leaf_shape,
                                    chunk_shape_with_color_as_tuple=chunk_shape_with_color_as_tuple,
                                    dtype=dtype,
                                    dataset=dataset,
                                    is_dataset_transposed=True)
                        #with Pool(16) as pool :
                        #    foo = list(tqdm.tqdm(pool.imap(f, tile_id_list), total=len(tile_id_list)))
                        # for tile_id in tqdm.tqdm(tile_id_list):
                        #     leaf_id_within_tile = tile_hash[tile_id]
                        #     f(tile_id, leaf_id_within_tile)
                        print('About to process %d tiles' % len(tile_id_list))
                        futures = client.map(two_arg_dump_single_tile_id, tile_id_list, leaf_ids_per_tile_list, retries=2)
                        progress(futures, notebook=False)  # need notebook=False when running in Spyder
                        wait(futures)  # just to make sure...
                        print('')
                        print('All Dask jobs have exited')
                        print('')
                        print('futures:')
                        print(futures)

                        #for tile_id in tile_id_list:
                        #    leaf_id_within_tile = tile_hash[tile_id]
                        #    this_future = client.submit(f, tile_id, leaf_id_within_tile)
                        #    fire_and_forget(this_future)
# end



# def sample_spherical(npoints, ndim=3):
#     vec = np.random.randn(ndim, npoints)
#     vec /= np.linalg.norm(vec, axis=0)
#     return vec



# def fixKinksinAnnotation():
#     input_folder = '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo'
#     swcfiles = [os.path.join(input_folder, fold, files) for fold in os.listdir(input_folder) if
#                 os.path.isdir(os.path.join(input_folder, fold)) for files in
#                 os.listdir(os.path.join(input_folder, fold)) if
#                 files.endswith("-carved.swc")]
#     swcfiles.sort()
#     for swc_file in swcfiles[1]:
#         path, filename = os.path.split(swc_file)
#         output_h5_file = os.path.join(path, filename.split('.')[0] + '.h5')
#         input_swc = swc_file
#         # output_h5_file = os.path.join(path,filename.split('.')[0][:-1]+'.h5')
#         # input_swc = os.path.join(path,filename.split('.')[0][:-1]+'.swc')
#         with h5py.File(output_h5_file, "r+") as f:
#             try:
#                 del f['reconstruction']
#             except Exception:
#                 pass
#             um, edges, R, offset, scale, header = util.readSWC(swcfile=input_swc, scale=1)
#             dset_swc = f.create_dataset("reconstruction", (um.shape[0], 7), dtype='f')
#             for iter, xyz_ in enumerate(um):
#                 xyz_ = np.ceil(xyz_ - np.sqrt(np.finfo(float).eps))
#                 dset_swc[iter, :] = np.array(
#                     [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])



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
        did_load_tile_hash = True
    except (OSError,FileNotFoundError) :
        did_load_tile_hash = False

    if not did_load_tile_hash :
        octpath_cover = np.unique(octpath, axis=0)
        #gridlist_cover = improc.oct2grid(octpath_cover)

        print('About to start dilation...')
        #octpath_dilated = octpath_cover.copy()
        desired_carve_out_half_diagonal_as_scalar = 512  # in x,y.  z will be different
        desired_carve_out_half_diagonal = desired_carve_out_half_diagonal_as_scalar * np.array([1.0, 1.0, spacing_um[0]/spacing_um[2]])
        #dilation_count = 8
        dilation_count = np.max( np.ceil(desired_carve_out_half_diagonal.astype(float) / leaf_shape.astype(float)) ).astype(int).item()
        # should be enough to get about a 512 vx cube around each swc centerpoint
        # (except 4x less in z, b/c axial rez is less)
        # t = time.time()
        # octpath_dilated_old = improc.dilateOct(octpath_cover, dilation_count)
        # elapsed = time.time() - t
        # print('Elapsed time for old method: %g s' % elapsed)
        t = time.time()
        octpath_dilated = improc.dilate_octree_chunk_set(octpath_cover, dilation_count)
        elapsed = time.time() - t
        print('Elapsed time for new dilation method: %g s' % elapsed)
        # if np.array_equal(octpath_dilated_old, octpath_dilated):
        #     print('The two methods agree on octpath dilation result!  Hooray!')
        # else:
        #     raise RuntimeError('The two methods do not agree on octpath dilation result')
        print('Done with dilation!')

        tile_hash = improc.chunklist(octpath_dilated, tile_level_count) #1..8
        os.makedirs(output_folder_name, exist_ok=True)
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
    dump_write(render_folder_name,
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



def main(argv):
    """ creates cropped volume and JW structure (for visualization) based on input render folder and swc file
        USAGE: 'navigator.py -i <data_folder> -s <swc_file> -o <output_folder>'
            -i <data_folder>: input data folder. Folders should follow octree format, e.g. <data_folder>/1/5/6
            -s <swc_file>: input swc_file or folder. for *swc files 7 column conventional reconstruction format.
            -o <output_folder>: folder to create h5 and JW files
            -h <number_of_level>: [OPTIONAL] sets how many chunks around trace will be used
            -j <output_octree>: [OPTIONAL] creates an octree formated folder at target location. "-j" without argument
                                will create target output at <output_folder>/JW location

        NOTES:
            oct in [1...8]
            grid in [0...(2**depth-1)]

            we keep mouselight data in <root>/<neuron-id>/consensus/<tag>_consensus.swc format, e.g.:
            /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-08-01/G-002/consensus/2018-08-01_G-002_consensus.swc
            it is suggested to copy all consensus files for that sample into a single folder manually or with a script than pass input folder with "-f" argument.
            For example:
            cd /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-08-01
            find . -name "*consensus*.swc" -exec cp {} /groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data/swc_recons/2018-08-01 \;


    """

    # @@TODO: multi swc data dump
    # @@TODO: fix octree dilation amount. make it user specified

    # data_fold='/nrs/mouselight/SAMPLES/2018-08-01-raw-rerender'
    # ## input_swc_file='/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo/2017-09-25_G-001_consensus/2017-09-25_G-001_consensus-proofed.swc'
    # input_swc_file='/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data/swc_recons/2018-08-01'
    # output_folder='/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo/2018-08-01'
    # octree_folder = os.path.join(output_folder, 'JW')

    additional_level_count = 3

    try:
        opts, args = getopt.getopt(argv, "i:s:o:j:f", ["data_fold=","input_swc_file=","output_folder=","octree_folder=","for"])
    except getopt.GetoptError:
        print('navigator.py -i <data_folder> -s <swc_file> -o <output_folder> -j <OPT:octree_folder>')
        sys.exit(2)


    do_use_simple_for_loop = False
    for opt, arg in opts:
        print('opt:', opt,'arg:', arg)
        if opt == '-h':
            print('navigator.py -i <data_folder> -s <swc_file> -o <output_folder>')
            sys.exit()
        elif opt in ("-i", "--data_fold"):
            print(arg)
            render_folder_name = arg
        elif opt in ("-s", "--input_swc_file"):
            input_swc_file = arg
        elif opt in ("-o", "--output_folder"):
            output_folder = arg
            octree_folder = os.path.join(output_folder,'JW')
        elif opt in ("-h", "--number_of_level"):
            additional_level_count = arg
        elif opt in ('-f', '--for'):
            do_use_simple_for_loop = True
        elif opt in ("-j", "--octree_folder"):
            try:
                octree_folder
            except NameError:
                print("Using output folder as JW folder")
                if octree_folder:
                    octree_folder = arg



    print('SWC FILE                  :', input_swc_file)
    print('DATA FOLDER               :', render_folder_name)
    print('OUTPUT FOLDER             :', output_folder)
    print('ADDITIONAL LEVEL COUNT    :', additional_level_count)
    print('do_use_simple_for_loop    :', do_use_simple_for_loop)
    #print('OCTREEFOLDER    :', octree_folder)


    rootfolder, swc_file_name = os.path.split(input_swc_file)
    #swc_name, _ = swc_name_w_ext.split(os.extsep)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    #output_swc_name = '{}-carved.swc'.format(swc_name)
    output_volume_file_name =  '{}-carved.n5'.format(swc_file_name)
    #JW_output_folder = os.path.join(output_folder,'JW')

    # if not os.path.exists(JW_output_folder):
    crop_from_render(render_folder_name, input_swc_file, output_folder, output_volume_file_name, do_use_simple_for_loop)

    # # shutil.rmtree(JW_output_folder)
    # if not os.path.exists(JW_output_folder):
    #     os.makedirs(JW_output_folder)
    #     os.chmod(JW_output_folder, 0o770)
    #
    # output_h5_file = os.path.join(output_folder, output_h5_name)
    # converter = util.Convert2JW(output_h5_file, JW_output_folder, number_of_oct_level=None)
    # converter.convert2JW()
    # converter.mergeJW(number_of_level=converter.number_of_oct_level)
    # converter.create_transform_file()

    print('DONE')
# end def



if __name__ == "__main__":
   main(sys.argv[1:])

