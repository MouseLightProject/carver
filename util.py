from functools import partial
#from multiprocessing import Pool
import numpy as np
import os
import h5py
import z5py
import tqdm
from dask_jobqueue import LSFCluster
from dask.distributed import Client, LocalCluster
#from dask.distributed import fire_and_forget
from dask.distributed import wait
from dask.distributed import progress
import getpass
#from collections import defaultdict
import improc

def readTransform(transformfile):
    # reads transform.txt file and parse it into a transform
    A = np.zeros((3,4))
    with open(transformfile, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            parts = text.split(':')
            num = np.float(parts[1].strip('\n'))
            if parts[0] == 'ox':
                A[0,3] = num
            elif parts[0] == 'oy':
                A[1, 3] = num
            elif parts[0] == 'oz':
                A[2,3] = num
            elif parts[0] == 'sx':
                A[0, 0] = num
            elif parts[0] == 'sy':
                A[1, 1] = num
            elif parts[0] == 'sz':
                A[2, 2] = num
            elif parts[0] == 'nl':
                # normalize diagonal with level
                np.fill_diagonal(A,A.diagonal()/(2**(num-1)))
    return A

def readParameterFile(parameterfile = ""):
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
    with open(parameterfile, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            parts = text.split('=')
            keyval = parts[0].strip()
            if keyval == 'const nlevels':
                params['nlevels'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const shape_leaf_px':
                params['leafSize'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const voxelsize_used_um':
                params['spacing'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const origin_nm':
                params['origin'] = 0.001 * np.array(eval(parts[1].strip('\n')),dtype=np.float)  # convert to um
            elif keyval == 'const nchannels':
                params['nchannels'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            else:
                it=0
    #A = np.zeros((3,4))
    #np.fill_diagonal(A, params['vixsize']*1000) #convert to nm
    #A[:,3] = params['origin']
    #params['A'] = A
    return params

def readSWCfolder(swcfolder, scale=1.0):
    swc_name_w_ext = os.listdir(swcfolder)
    nm, edges, R, offset, header, filenames = ([] for i in range(6))
    curr_len = 0
    for iswc_name in swc_name_w_ext:
        nm_, edges_, R_, offset_, scale_, header_, junk = readSWC(swcfile=os.path.join(swcfolder, iswc_name), scale=scale)

        nm.append(nm_)
        edges.append(edges_)
        R.append(R_)
        offset.append(offset_)
        header.append(header_)
        filenames.append(iswc_name)

    return nm, edges, R, filenames, header

def appendSWCfolder(swcfolder, scale=1.0):
    swc_name_w_ext = os.listdir(swcfolder)

    nm, edges, R, offset, header, filenames = ([] for i in range(6))

    curr_len = 0
    for iswc_name in swc_name_w_ext:
        nm_, edges_, R_, offset_, scale_, header_, junk = readSWC(swcfile=os.path.join(swcfolder, iswc_name), scale=scale)
        edges_[0,1] = edges_[0,0] # swc convention sets root to -1
        edges_ += curr_len # append to previous recon
        curr_len += nm_.shape[0]

        nm.append(nm_)
        edges.append(edges_)
        R.append(R_)
        offset.append(offset_)
        header.append(header_)
        filenames.append(iswc_name)

    # concatenate
    nm_ = np.vstack(nm)
    edges_ = np.vstack(edges)
    R_ = np.hstack(R)
    return nm_, edges_, R_


def collect_nodes_from_tracing_complete_folder(folder_name):
    # Get names of files and folders
    folder_entries = os.listdir(folder_name)

    # only want folder entries of the form G-%03d
    neuron_names = [entry for entry in folder_entries if entry.startswith('G-') and len(entry)==5]

    #xyzs = ([] for i in range(1))
    xyzs = []
    for neuron_name in neuron_names:
        #xyz_this, edges_, R_, offset_, scale_, header_, junk = readSWC(swcfile=os.path.join(folder_name, iswc_name))
        xyzs_this = readSWC(os.path.join(folder_name, neuron_name, 'consensus.swc'))[0]
        xyzs.append(xyzs_this)
        xyzs_this = readSWC(os.path.join(folder_name, neuron_name, 'dendrite.swc'))[0]
        xyzs.append(xyzs_this)

    # Make the list of lists into a single list
    xyzs_stacked = np.vstack(xyzs)
    
    return xyzs_stacked


# ORIGINAL_SOURCE Janelia Workstation Large Volume Viewer
# OFFSET 66310.961575 46976.514329 18608.718278
# COLOR 1.000000,0.200000,0.200000
def readSWC(swcfile, scale=1.0):
    swcline=[]
    offset = np.zeros((1,3))
    offsetkey = 'OFFSET'
    header = []
    with open(swcfile, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            if text[0]=='#':
                header.append(text)
                # check offset
                if text[2:len(offsetkey)+2]==offsetkey:
                    offset = np.array(text[len(offsetkey) + 3:-1].split(), dtype=np.float).reshape(1,3)
                else:
                    continue #skip
            else:
                parts = text.split()
                swcline.append(parts)
    lines = np.array(swcline, dtype=float).reshape(-1, 7)
    edges = lines[:,(0,6)]
    R = lines[:,5]
    xyz = lines[:,2:5]
    xyz = xyz + offset
    xyz = xyz/scale
    structure_identifier = lines[:,1]
    return (xyz,edges,R,offset,scale,header,structure_identifier)

def upsampleSWC(xyz, edges, sp):
    if xyz.shape[0]==1:
        return xyz
    xyzup = []
    for i, j in np.asarray(edges - 1, np.int64):
        if j < 0:
            continue
        else:
            st = xyz[i, :][None, :]
            ed = xyz[j, :][None, :]
            el = ed - st
            enel = el / np.linalg.norm(el)
            numiter = np.ceil(np.linalg.norm(el) / sp)
            xyzup.append(np.arange(0, numiter).reshape(-1, 1) * enel * sp + st)

    return(np.concatenate(xyzup))


def um2pix(xyz_um, origin_um, spacing_um):
    # Applies transform to convert um into pix location
    return( np.round((xyz_um-origin_um)/spacing_um) )


def pix2um(xyz_voxels, origin_um, spacing_um):
    return( spacing_um*xyz_voxels + origin_um )


# def dump_single_tile_id(tile_id, 
#                         leaf_ids_within_tile, 
#                         rendered_folder_path, 
#                         tile_shape, 
#                         leaf_shape, 
#                         chunk_shape_with_color_as_tuple, 
#                         dtype, 
#                         dataset, 
#                         is_dataset_transposed):
#     tilename = '/'.join(a for a in tile_id)
#     tilepath = os.path.join(rendered_folder_path, tilename)

#     tile_octree_path = np.array(list(tile_id), dtype=int)
#     tile_ijk_in_tile_grid = np.ndarray.flatten(improc.oct2grid(tile_octree_path.reshape(1, len(tile_octree_path))))
#     tile_origin_ijk = tile_ijk_in_tile_grid * tile_shape
#     tile_end_ijk = tile_origin_ijk + tile_shape
#     if os.path.isdir(tilepath):
#         tile_ijk_in_tile_grid_as_tuple = tuple(tile_ijk_in_tile_grid)
#         tile_ijk_in_tile_grid_as_tuple_with_color = tile_ijk_in_tile_grid_as_tuple + (0,)  # dataset includes color channels
#         if is_dataset_transposed:
#             tile_ijk_in_tile_grid_as_tuple_with_color_maybe_flipped = tuple(reversed(tile_ijk_in_tile_grid_as_tuple_with_color)) 
#         else:
#             tile_ijk_in_tile_grid_as_tuple_with_color_maybe_flipped = tile_ijk_in_tile_grid_as_tuple_with_color
#         does_chunk_exist = dataset.chunk_exists(tile_ijk_in_tile_grid_as_tuple_with_color_maybe_flipped)
#         if not does_chunk_exist:
#             im = improc.loadTiles(tilepath)
#             # relativeDepth = leaf_level_count - tile_level_count
#             output_tile_stack = np.zeros(chunk_shape_with_color_as_tuple, dtype=dtype)

#             # patches in idTiled
#             for leaf_octree_path_within_tile_as_string in leaf_ids_within_tile:
#                 leaf_octree_path_within_tile = np.array(list(leaf_octree_path_within_tile_as_string), dtype=int)
#                 leaf_ijk_in_leaf_grid_within_tile = improc.oct2grid(leaf_octree_path_within_tile.reshape(1, len(leaf_octree_path_within_tile)))  # in 0 base

#                 start = np.ndarray.flatten(leaf_ijk_in_leaf_grid_within_tile * leaf_shape)
#                 end = np.ndarray.flatten(start + leaf_shape)
#                 leaf_stack = im[start[0]:end[0], start[1]:end[1], start[2]:end[2], :]
#                 output_tile_stack[start[0]:end[0], start[1]:end[1], start[2]:end[2], :] = leaf_stack

#             if is_dataset_transposed :
#                 dataset[:, tile_origin_ijk[2]:tile_end_ijk[2], tile_origin_ijk[1]:tile_end_ijk[1], tile_origin_ijk[0]:tile_end_ijk[0]] = \
#                     np.transpose(output_tile_stack)                
#             else:
#                 dataset[tile_origin_ijk[0]:tile_end_ijk[0], tile_origin_ijk[1]:tile_end_ijk[1], tile_origin_ijk[2]:tile_end_ijk[2], :] = \
#                     output_tile_stack


# def dump_write(render_folder_name,
#                full_volume_shape,
#                dtype,
#                color_channel_count,
#                output_file_name,
#                tile_hash,
#                leaf_level_count,
#                tile_level_count,
#                compression_method,
#                compression_options,
#                output_file_type,
#                do_use_simple_for_loop=False):
#     # dumps volumetric data into h5/n5/zarr
#     #self.inputLoc = inputloc

#     tile_shape = (full_volume_shape / (2**tile_level_count)).astype(int)
#     leaf_shape = (full_volume_shape / (2**leaf_level_count)).astype(int)

#     # check if dataset name is provided
#     splitted_name = output_file_name.split(':')
#     if  len(splitted_name) == 1:
#         output_file_name =  splitted_name[0]
#         dataset_name =  "volume"
#     elif len(splitted_name) ==2:
#         output_file_name =  splitted_name[0]
#         dataset_name =  splitted_name[1]
#     else:
#         raise ValueError('output file name has more than one ":"', output_file_name)
#     #self.setting = setting
#     #self.tilelist = tilelist
#     tile_id_list = list(tile_hash.keys())
#     leaf_ids_per_tile_list = list(tile_hash.values())

#     # # Unpack the settings
#     # volSize = tuple(map(int,setting['volSize']))
#     # tileSize = setting['tileSize']
#     # #volReference = setting['volReference']
#     # depthFull = setting['depthFull']
#     # depthBase = setting['depthBase']
#     # leafSize = setting['leaf_shape']
#     # dtype = setting['dtype']
#     # chunkSize = tuple(map(int,setting['chunkSize']))
#     # compression_method = setting['compression']
#     # comp_opts = setting['compression_opts']
#     chunk_shape = tile_shape
#     full_volume_shape_including_color_channel = np.append(full_volume_shape, color_channel_count)  # append color channel
#     chunk_shape_including_color_channel = np.append(chunk_shape, color_channel_count)
#     full_volume_shape_with_color_channels_as_tuple = tuple(map(int, full_volume_shape_including_color_channel))
#     chunk_shape_with_color_as_tuple = tuple(map(int, chunk_shape_including_color_channel))

#     if output_file_type=='h5':
#         # write into h5
#         with h5py.File(output_file_name, "w") as f:
#             # dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0], 7), dtype='f')
#             # for iter, xyz_ in enumerate(xyz_shifted):
#             #     dset_swc[iter, :] = np.array(
#             #         [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])
#             dataset = f.create_dataset(dataset_name,
#                                        full_volume_shape_with_color_channels_as_tuple,
#                                        dtype=dtype,
#                                        chunks=chunk_shape_with_color_as_tuple,
#                                        compression=compression_method,
#                                        compression_opts=compression_options)


#             # crop chuncks from a tile read in tilelist
#             for iter, tile_id in enumerate(tile_id_list):
#                 print('{} : {} out of {}'.format(tile_id, iter+1, len(tile_id_list)))
#                 leaf_id_within_tile = tile_hash[tile_id]
#                 dump_single_tile_id(tile_id,
#                                     leaf_id_within_tile,
#                                     render_folder_name,
#                                     tile_shape,
#                                     leaf_shape,
#                                     chunk_shape_with_color_as_tuple,
#                                     dtype,
#                                     dataset, 
#                                     is_dataset_transposed=False)
#     elif output_file_type=='n5' or output_file_type=='zarr':
#         # write into z5 or n5
#         if do_use_simple_for_loop:
#             use_zarr_format = (output_file_type == 'zarr')
#             with z5py.File(output_file_name, 'a', use_zarr_format=use_zarr_format) as f:
#                 # require_dataset seems to choke on the compression_options {level: 9}, so this is a workaround
#                 g = f.require_group('/')
#                 try:
#                     dataset = g[dataset_name]
#                 except KeyError:
#                     dataset = f.create_dataset(dataset_name,
#                                                shape=tuple(reversed(full_volume_shape_with_color_channels_as_tuple)),
#                                                dtype=dtype,
#                                                chunks=tuple(reversed(chunk_shape_with_color_as_tuple)),
#                                                compression=compression_method,
#                                                **compression_options)                    
#                 for tile_id in tqdm.tqdm(tile_id_list):
#                     leaf_ids_within_tile = tile_hash[tile_id]
#                     dump_single_tile_id(tile_id,
#                                         leaf_ids_within_tile,
#                                         render_folder_name,
#                                         tile_shape,
#                                         leaf_shape,
#                                         chunk_shape_with_color_as_tuple,
#                                         dtype,
#                                         dataset,
#                                         is_dataset_transposed=True)
#         else:
#             username = getpass.getuser()
#             scratch_folder_path = '/scratch/%s' % username
#             with LSFCluster(cores=1, memory='15 GB', local_dir=scratch_folder_path, projectstr='mouselight', queue='normal', extralist='-o /dev/null -e /dev/null') as cluster:
#                 #cluster.adapt(minimum=1, maximum=1000)
#                 #cluster = LocalCluster(n_workers=4, threads_per_worker=1)
#                 cluster.scale(200)
#                 with Client(cluster) as client:
#                     use_zarr_format = (output_file_type=='zarr')
#                     with z5py.File(output_file_name, 'a', use_zarr_format=use_zarr_format) as f:
#                         # require_dataset seems to choke on the compression_options {level: 9}, so this is a workaround
#                         g = f.require_group('/')
#                         try:
#                             dataset = g[dataset_name]
#                         except KeyError:                        
#                             dataset = f.create_dataset(dataset_name,
#                                                        shape=tuple(reversed(full_volume_shape_with_color_channels_as_tuple)),
#                                                        dtype=dtype,
#                                                        chunks=tuple(reversed(chunk_shape_with_color_as_tuple)),
#                                                        compression=compression_method,
#                                                        **compression_options)
#                         two_arg_dump_single_tile_id = \
#                             partial(dump_single_tile_id,
#                                     rendered_folder_path=render_folder_name,
#                                     tile_shape=tile_shape,
#                                     leaf_shape=leaf_shape,
#                                     chunk_shape_with_color_as_tuple=chunk_shape_with_color_as_tuple,
#                                     dtype=dtype,
#                                     dataset=dataset,
#                                     is_dataset_transposed=True)
#                         #with Pool(16) as pool :
#                         #    foo = list(tqdm.tqdm(pool.imap(f, tile_id_list), total=len(tile_id_list)))
#                         # for tile_id in tqdm.tqdm(tile_id_list):
#                         #     leaf_id_within_tile = tile_hash[tile_id]
#                         #     f(tile_id, leaf_id_within_tile)
#                         print('About to process %d tiles' % len(tile_id_list))
#                         futures = client.map(two_arg_dump_single_tile_id, tile_id_list, leaf_ids_per_tile_list, retries=2)
#                         progress(futures, notebook=False)  # need notebook=False when running in Spyder
#                         wait(futures)  # just to make sure...
#                         print('')
#                         print('All Dask jobs have exited')
#                         print('')
#                         print('futures:')
#                         print(futures)

#                         #for tile_id in tile_id_list:
#                         #    leaf_id_within_tile = tile_hash[tile_id]
#                         #    this_future = client.submit(f, tile_id, leaf_id_within_tile)
#                         #    fire_and_forget(this_future)


