import util
import numpy as np
import matplotlib.pyplot as plt
import improc
import os
from scipy.spatial import distance as dist
import h5py

from mpl_toolkits.mplot3d import Axes3D

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def crop_from_render():
    # decomp
    # tif2mj2(inputtif,mj2file)
    test = 1

if __name__ == '__main__':

    # oct in [1...8]
    # grid in [0...(2**depth-1)]

    epsball = 500

    import importlib

    importlib.reload(util)
    importlib.reload(improc)
    datafold = '/nrs/mouselight/SAMPLES/2017-09-11-2'

    params = util.readParameterFile(parameterfile=datafold+"/calculated_parameters.jl")
    um, edges, R, offset, scale = util.readSWC(swcfile='./2017-09-11_G-001_CA_base.swc',scale=1/1000)
    # to fix the bug in Janelia Workstation
    um = um + params['vixsize']/scale/2
    xyz = util.um2pix(um,params['A']).T
    # upsample xyz to
    sp = 5

    xyzup = util.upsampleSWC(xyz, edges, sp)
    if False:
        octpath, xres = improc.xyz2oct(xyzup,params)
    else:
        depthextend = 3
        params_p1=params.copy()
        params_p1["nlevels"]=params_p1["nlevels"]+depthextend
        params_p1["leafshape"]=params_p1["leafshape"]/(2**depthextend)
        octpath, xres = improc.xyz2oct(xyzup,params_p1)

    depthBase = params["nlevels"].astype(int)
    depthFull = params_p1["nlevels"].astype(int)
    tileSize = params["leafshape"].astype(int)
    leafSize = params_p1["leafshape"].astype(int)

    octpath_cover = np.unique(octpath, axis=0)
    gridlist_cover = improc.oct2grid(octpath_cover)
    octpath_dilated,gridlist_dilated = improc.dilateOct(octpath_cover)

    tilelist = improc.chunklist(octpath_dilated,depthBase) #1..8
    tileids = list(tilelist.keys())

    gridReference = np.min(gridlist_dilated, axis=0)
    volReference = gridReference*leafSize
    gridSize = np.max(gridlist_dilated, axis=0) - gridReference +1
    outVolumeSize = np.append(gridSize*leafSize,2) #append color channel
    chunksize = np.append(leafSize,2)

    # save into cropped swc
    xyz_shifted = xyz-volReference
    with open('test.swc','w') as fswc:
        for iter,txt in enumerate(xyz_shifted):
            fswc.write('{:d} {} {:.2f} {:.2f} {:.2f} {:.2f} {:d}\n'.format(edges[iter,0].__int__(),1,txt[0],txt[1],txt[2],1,edges[iter,1].__int__()))

    # write into h5
    with h5py.File("test2.hdf5", "w") as f:
        dset = f.create_dataset("volume", outVolumeSize, dtype='uint16', chunks=tuple(chunksize), compression="gzip", compression_opts=9)
        # crop chuncks from a tile read in tilelist
        for iter,idTile in enumerate(tileids):
            print('{} : {} out of {}'.format(idTile, iter, len(tileids)))
            tilename = '/'.join(a for a in idTile)
            tilepath = datafold+'/'+tilename

            ijkTile = np.array(list(idTile), dtype=int)
            xyzTile = improc.oct2grid(ijkTile.reshape(1, len(ijkTile)))
            locTile = xyzTile * tileSize
            locShift = np.asarray(locTile - volReference,dtype=int).flatten()
            im = improc.loadTiles(tilepath)
            relativeDepth = depthFull - depthBase

            # patches in idTiled
            for patch in tilelist[idTile]:
                ijk = np.array(list(patch),dtype=int)
                xyz = improc.oct2grid(ijk.reshape(1, len(ijk))) # in 0 base

                start = np.ndarray.flatten(xyz*leafSize)
                end = np.ndarray.flatten(start + leafSize)
                print(start,end)
                imBatch = im[start[0]:end[0],start[1]:end[1],start[2]:end[2],:]

                # # viz patch
                # fig = plt.figure()
                # plt.imshow(np.max(imBatch[:,:,:,0],axis=2))

                start = start + locShift
                end = end + locShift
                dset[start[0]:end[0],start[1]:end[1],start[2]:end[2],:] = imBatch

        dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0],7), dtype='f')
        for iter, xyz_ in enumerate(xyz_shifted):
            dset_swc[iter,:] = np.array([edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1, edges[iter, 1].__int__()])
