import numpy as np
import itertools
import os
from skimage import io
import h5py
from skimage.transform import resize
import warnings
import z5py

from collections import defaultdict
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
        nm_, edges_, R_, offset_, scale_, header_ = readSWC(swcfile=os.path.join(swcfolder, iswc_name), scale=scale)

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
        nm_, edges_, R_, offset_, scale_, header_ = readSWC(swcfile=os.path.join(swcfolder, iswc_name), scale=scale)
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
    return (xyz,edges,R,offset,scale,header)

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
    # applies transform to convert um into pix location
    # um_ = np.concatenate((um, np.ones((um.shape[0], 1))), axis=1)
    # return(np.dot(np.linalg.pinv(A),um.T))
    #return(np.dot(np.diag(1 / np.diagonal(A[:3, :3])), (xyz_um - A[:, 3]).T))
    return( (xyz_um-origin_um)/spacing_um )


def pix2um(xyz_voxels, origin_um, spacing_um):
    #return(np.dot(A, xyz_voxels))
    return( spacing_um*xyz_voxels + origin_um )

# def pix2oct(xyz,dims,depth):
#     # for a given xyz, box size and depth, returns the location int the patch and patch path
#     res = dims/depth
#     ijk = np.floor(xyz/res)
#     # convert ijk to
#     return 0
#
#
# def um2oct(xyz,dims,transform ):
#     # for a given um, transform and image size, returns the patch location
#     return 0
#
#
# def traverseOct():
#     # lets you to traverse octree
#     return 0


# class dumper(object):
#     # dumps volumetric data into h5/zarr
#     def __init__(self, inputloc, outputFile, setting, tilelist=None):
#         self.inputLoc = inputloc
#         # check if dataset name is provided
#         splitted_name = outputFile.split(':')
#         if  len(splitted_name) == 1:
#             self.outputFile =  splitted_name[0]
#             self.datasetName =  "volume"
#         elif len(splitted_name) ==2:
#             self.outputFile =  splitted_name[0]
#             self.datasetName =  splitted_name[1]
#         else:
#             raise ValueError('output file name has more than one ":"', outputFile)
#         self.setting = setting
#         self.tilelist = tilelist
#         if tilelist:
#             self.tileids = list(tilelist.keys())
#
#     def write(self):
#         if self.setting['type'] is 'h5':
#             # write into h5
#             tileids = self.tileids
#             inputLoc = self.inputLoc
#             outputFile = self.outputFile
#             tilelist = self.tilelist
#             setting = self.setting
#             volSize = setting['volSize']
#             tileSize = setting['tileSize']
#             volReference = setting['volReference']
#             depthFull = setting['depthFull']
#             depthBase = setting['depthBase']
#             leafSize = setting['leafSize']
#
#             with h5py.File(outputFile, "a") as f:
#                 # dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0], 7), dtype='f')
#                 # for iter, xyz_ in enumerate(xyz_shifted):
#                 #     dset_swc[iter, :] = np.array(
#                 #         [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])
#                 dset = f.create_dataset(self.datasetName, volSize, dtype=setting['dtype'], chunks=setting['chunkSize'],
#                                         compression=setting['compression'], compression_opts=setting['compression_opts'])
#                 # crop chuncks from a tile read in tilelist
#                 for iter, idTile in enumerate(tileids):
#                     print('{} : {} out of {}'.format(idTile, iter, len(tileids)))
#                     tilename = '/'.join(a for a in idTile)
#                     tilepath = os.path.join(inputLoc, tilename)
#
#                     ijkTile = np.array(list(idTile), dtype=int)
#                     xyzTile = improc.oct2grid(ijkTile.reshape(1, len(ijkTile)))
#                     locTile = xyzTile * tileSize
#                     locShift = np.asarray(locTile - volReference, dtype=int).flatten()
#                     if os.path.isdir(tilepath):
#
#                         im = improc.loadTiles(tilepath)
#                         relativeDepth = depthFull - depthBase
#
#                         # patches in idTiled
#                         for patch in tilelist[idTile]:
#                             ijk = np.array(list(patch), dtype=int)
#                             xyz = improc.oct2grid(ijk.reshape(1, len(ijk)))  # in 0 base
#
#                             start = np.ndarray.flatten(xyz * leafSize)
#                             end = np.ndarray.flatten(start + leafSize)
#                             # print(start,end)
#                             imBatch = im[start[0]:end[0], start[1]:end[1], start[2]:end[2], :]
#
#                             start = start + locShift
#                             end = end + locShift
#                             dset[start[0]:end[0], start[1]:end[1], start[2]:end[2], :] = imBatch


def dump_write(inputLoc, outputFile, setting, tilelist):
    # dumps volumetric data into h5/zarr
    #self.inputLoc = inputloc
    # check if dataset name is provided
    splitted_name = outputFile.split(':')
    if  len(splitted_name) == 1:
        outputFile =  splitted_name[0]
        datasetName =  "volume"
    elif len(splitted_name) ==2:
        outputFile =  splitted_name[0]
        datasetName =  splitted_name[1]
    else:
        raise ValueError('output file name has more than one ":"', outputFile)
    #self.setting = setting
    #self.tilelist = tilelist
    tileids = list(tilelist.keys())

    # Unpack the settings
    volSize = tuple(map(int,setting['volSize']))
    tileSize = setting['tileSize']
    volReference = setting['volReference']
    depthFull = setting['depthFull']
    depthBase = setting['depthBase']
    leafSize = setting['leaf_shape']
    dtype = setting['dtype']
    chunkSize = tuple(map(int,setting['chunkSize']))
    compression_method = setting['compression']
    comp_opts = setting['compression_opts']

    if setting['type'] is 'h5':
        # write into h5
        with h5py.File(outputFile, "w") as f:
            # dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0], 7), dtype='f')
            # for iter, xyz_ in enumerate(xyz_shifted):
            #     dset_swc[iter, :] = np.array(
            #         [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])
            dset = f.create_dataset(datasetName,
                                    volSize,
                                    dtype=dtype,
                                    chunks=chunkSize,
                                    compression=compression_method,
                                    compression_opts=comp_opts)
            # crop chuncks from a tile read in tilelist
            for iter, idTile in enumerate(tileids):
                print('{} : {} out of {}'.format(idTile, iter+1, len(tileids)))
                tilename = '/'.join(a for a in idTile)
                tilepath = os.path.join(inputLoc, tilename)

                ijkTile = np.array(list(idTile), dtype=int)
                xyzTile = improc.oct2grid(ijkTile.reshape(1, len(ijkTile)))
                locTile = xyzTile * tileSize
                locShift = np.asarray(locTile - volReference, dtype=int).flatten()
                if os.path.isdir(tilepath):

                    im = improc.loadTiles(tilepath)
                    relativeDepth = depthFull - depthBase

                    # patches in idTiled
                    for patch in tilelist[idTile]:
                        ijk = np.array(list(patch), dtype=int)
                        xyz = improc.oct2grid(ijk.reshape(1, len(ijk)))  # in 0 base

                        start = np.ndarray.flatten(xyz * leafSize)
                        end = np.ndarray.flatten(start + leafSize)
                        # print(start,end)
                        imBatch = im[start[0]:end[0], start[1]:end[1], start[2]:end[2], :]

                        start = start + locShift
                        end = end + locShift
                        dset[start[0]:end[0], start[1]:end[1], start[2]:end[2], :] = imBatch
    elif setting['type'] is 'n5' or setting['type'] is 'zarr':
        # write into z5 or n5
        use_zarr_format = (setting['type']=='zarr')
        with z5py.File(outputFile, "w", use_zarr_format=use_zarr_format) as f:
            dset = f.create_dataset(datasetName,
                                    shape=volSize,
                                    dtype=dtype,
                                    chunks=chunkSize,
                                    compression=compression_method,
                                    **comp_opts)
            # crop chuncks from a tile read in tilelist
            for iter, idTile in enumerate(tileids):
                print('{} : {} out of {}'.format(idTile, iter+1, len(tileids)))
                tilename = '/'.join(a for a in idTile)
                tilepath = os.path.join(inputLoc, tilename)

                ijkTile = np.array(list(idTile), dtype=int)
                xyzTile = improc.oct2grid(ijkTile.reshape(1, len(ijkTile)))
                locTile = xyzTile * tileSize
                locShift = np.asarray(locTile - volReference, dtype=int).flatten()
                if os.path.isdir(tilepath):

                    im = improc.loadTiles(tilepath)
                    relativeDepth = depthFull - depthBase

                    # patches in idTiled
                    for patch in tilelist[idTile]:
                        ijk = np.array(list(patch), dtype=int)
                        xyz = improc.oct2grid(ijk.reshape(1, len(ijk)))  # in 0 base

                        start = np.ndarray.flatten(xyz * leafSize)
                        end = np.ndarray.flatten(start + leafSize)
                        # print(start,end)
                        imBatch = im[start[0]:end[0], start[1]:end[1], start[2]:end[2], :]

                        start = start + locShift
                        end = end + locShift
                        dset[start[0]:end[0], start[1]:end[1], start[2]:end[2], :] = imBatch



class Convert2JW(object):
    def __init__(self,h5file,experiment_folder,number_of_oct_level=None):
        with h5py.File(h5file, "r") as f:
            volume = f["volume"]
            self.h5_dims= np.array(volume.shape)
            self.h5_chunk_size = np.array(volume.chunks)
            if not number_of_oct_level:
                # estimate leaf size & depth
                # use multiple of chunk size for leaf
                self.target_leaf_size = self.h5_chunk_size[:3] * 8  # set target_leaf_size to a multiple of chunk size for efficiency
                depths = np.arange(2, 7)[:,None]
                self.output_dims = 2**depths[np.where(np.all(2**depths*self.target_leaf_size[None,:]>=self.h5_dims[:3],axis=1))[0]].flatten()[0]*self.target_leaf_size
                # depths = np.arange(2, 10)[:,None]
                # self.output_dims = 2**depths[np.where(np.all(2 **depths*self.h5_chunk_size[:3][None,:]>self.h5_dims[:3],axis=1))[0]].flatten()[0]*self.h5_chunk_size[:3]
                self.number_of_oct_level = np.log2(self.output_dims[0]/self.target_leaf_size[0]).__int__()
            else:
                self.output_dims = self.h5_dims
                self.number_of_oct_level = number_of_oct_level
                self.target_leaf_size = np.asarray(np.ceil(np.array(self.output_dims[:3]) / 2 ** self.number_of_level),
                                                   np.int)  # need to iter over leafs
        self.h5file = h5file
        self.experiment_folder = experiment_folder

    def __str__(self):
        return "{}\n{}\n{}\n{}\n{}".format(self.output_dims,self.number_of_level,self.target_leaf_size,self.h5file,self.experiment_folder)

    def convert2JW(self):
        h5file = self.h5file
        number_of_oct_level = self.number_of_oct_level
        experiment_folder = self.experiment_folder
        target_leaf_size = self.target_leaf_size
        with h5py.File(h5file, "r") as f:
            volume = f["volume"]
            output_dims = volume.shape
            bit_multiplication_array = 2**np.arange(3)
            #target_leaf_size = np.asarray(np.ceil(np.array(output_dims[:3]) / 2 ** number_of_oct_level), np.int)
            padded_size = target_leaf_size*2**number_of_oct_level
            range_values = [np.asarray(np.arange(0, padded_size[ii], target_leaf_size[ii]),dtype=np.int).tolist() for ii in range(3)]

            for ix,ref in enumerate(list(itertools.product(*range_values))):
                bb_end = np.asarray(np.min((ref+target_leaf_size,np.array(output_dims[:3])),axis=0),dtype=np.int)
                patch_ = volume[ref[0]:bb_end[0], ref[1]:bb_end[1], ref[2]:bb_end[2], :]
                if ~np.any(patch_):
                    continue

                # '''if patch size is smaller than full volume size pad zeros'''
                if np.any(bb_end-np.array(ref) < target_leaf_size):
                    # pad
                    patch = np.zeros(np.append(target_leaf_size,2),dtype=np.uint16)
                    patch[:patch_.shape[0], :patch_.shape[1], :patch_.shape[2],:] = patch_
                else:
                    patch = patch_

                folder_inds = np.array(np.unravel_index(ix, ([2**number_of_oct_level for ii in range(3)])))
                folder_inds = folder_inds + 1

                patch_folder_path = []
                for im in np.arange(number_of_oct_level,0,-1):
                    bits = folder_inds>2**(im-1)
                    # bit 2 num
                    patch_folder_path.append(1+np.sum(bits*bit_multiplication_array))
                    folder_inds = folder_inds-2**(im-1)*bits

                # create folder
                outfolder = os.path.join(experiment_folder,'/'.join(str(pp) for pp in patch_folder_path))
                # if ~np.any(patch):
                #     continue

                if not os.path.exists(outfolder):
                    os.makedirs(outfolder)
                print(outfolder)
                for ichannel in range(2):
                    outfile = os.path.join(outfolder,'default.'+str(ichannel)+'.tif')
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        io.imsave(outfile, np.swapaxes(patch[:, :, :, ichannel], 2, 0))

    def create_transform_file(self):
        experiment_folder = self.experiment_folder
        number_of_oct_level = self.number_of_oct_level
        # create transform file
        transform_file = os.path.join(experiment_folder,'transform.txt')
        with open(transform_file, 'w') as ft:
            ft.write('ox: {:04.0f}\n'.format(0))
            ft.write('oy: {:04.0f}\n'.format(0))
            ft.write('oz: {:04.0f}\n'.format(0))
            ft.write('sx: {:.0f}\n'.format(2**number_of_oct_level*1000))
            ft.write('sy: {:.0f}\n'.format(2**number_of_oct_level*1000))
            ft.write('sz: {:.0f}\n'.format(2**number_of_oct_level*1000))
            ft.write('nl: {:.0f}\n'.format(number_of_oct_level+1))

    def create_yml_file(self):
        yml_file = os.path.join(self.experiment_folder,'tilebase.cache.yml')

    def mergeJW(self,number_of_level=3):
        # reads an octant, down samples it and save a tif file
        # create all paths from depth-1 to 0
        experiment_folder = self.experiment_folder
        leaf_size = self.target_leaf_size
        values = ['{}/'.format(str(ii + 1)) for ii in range(8)]
        for current_number_of_level in np.arange(number_of_level-1, -1, -1):
            for iter_current_folder in list(itertools.product(values, repeat=current_number_of_level)):
                my_lst_str = ''.join(map(str, iter_current_folder))
                current_folder = os.path.join(experiment_folder, my_lst_str)
                if os.path.exists(current_folder):
                    print(current_folder)
                    self.__iter_octant__(current_folder,leaf_size)

    def __iter_octant__(self,current_folder,leaf_shape):
        for ichannel in range(2):
            im_channel = []
            for ioct in range(8):
                current_path = os.path.join(current_folder,str(ioct+1))
                current_file = current_path+'/default.{}.tif'.format(ichannel)
                if os.path.exists(current_file):
                    im_batch = np.swapaxes(io.imread(current_file), 2, 0)
                else:
                    im_batch = np.zeros(leaf_shape)
                im_channel.append(im_batch)

            rt1 = np.concatenate(im_channel[0:2], axis=0)
            rt2 = np.concatenate(im_channel[2:4], axis=0)
            rt3 = np.concatenate(im_channel[4:6], axis=0)
            rt4 = np.concatenate(im_channel[6:8], axis=0)
            rt5 = np.concatenate((rt1, rt2), axis=1)
            rt6 = np.concatenate((rt3, rt4), axis=1)
            merged_Im = np.concatenate((rt5, rt6), axis=2)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # down sample image by 2
                down_image = resize(merged_Im, leaf_shape, preserve_range=True)
                io.imsave(current_folder+'/default.{}.tif'.format(ichannel),
                          np.asarray(np.swapaxes(down_image, 2, 0), np.uint16))

