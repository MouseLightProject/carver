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

    params = util.readParameterFile(parameterfile="/nrs/mouselight/SAMPLES/2017-06-10/calculated_parameters.jl")
    um, edges, R, offset, scale = util.readSWC(scale=1/1000)
    # to fix the bug in Janelia Workstation
    um = um + params['vixsize']/scale/2
    xyz = util.um2pix(um,params['A']).T
    # upsample xyz to
    sp = 5

    xyzup = util.upsampleSWC(xyz, edges, sp)
    if False:
        octpath, xres = improc.xyz2oct(xyzup,params)
    else:
        depthextend = 1
        params_p1=params.copy()
        params_p1["nlevels"]=params_p1["nlevels"]+depthextend
        params_p1["leafshape"]=params_p1["leafshape"]/(2**depthextend)
        octpath, xres = improc.xyz2oct(xyzup,params_p1)

    depthBase = params["nlevels"].__int__()
    depthFull = params_p1["nlevels"].__int__()
    tileSize = params["leafshape"].astype(int)
    leafSize = params_p1["leafshape"].astype(int)
    octpath_cover = np.unique(octpath, axis=0)


    gridlist_cover = improc.oct2grid(octpath_cover)
    octpath_dilated,gridlist_dilated = improc.dilateOct(octpath_cover)

    tilelist = improc.chunklist(octpath_dilated,depthBase.__int__())
    tileids = list(tilelist.keys())

    gridReference = np.min(gridlist_dilated, axis=0)
    volReference = gridReference*leafSize
    gridSize = np.max(gridlist_dilated, axis=0) - gridReference +1
    outVolumeSize = np.append(gridSize*leafSize,2) #append color channel
    chunksize = np.append(leafSize,2)

    f = h5py.File("test.hdf5", "w")
    dset = f.create_dataset("volume", outVolumeSize, chunks=tuple(chunksize), compression="gzip", compression_opts=9)
    # f.__delitem__("volume")

    datafold = '/Volumes/mouselight/SAMPLES/2017-06-10'
    for idTile in tileids2:
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
            xyz = improc.oct2grid(ijk.reshape(1, len(ijk)))

            start = np.ndarray.flatten(xyz*leafSize)
            end = np.ndarray.flatten(start + leafSize)
            imBatch = im[start[0]:end[0],start[1]:end[1],start[2]:end[2],:]
            print(start,end)

            start = start + locShift
            end = end + locShift
            dset[start[0]:end[0],start[1]:end[1],start[2]:end[2],:] = imBatch


    # # for every batch find relative shift in coordinate
    # # reference shift
    # octpath_dilated_base = octpath_dilated.copy()
    # octpath_dilated_base[:,depthBase::]=1
    # gridlist_dilated_base = improc.oct2grid(octpath_dilated_base)
    # # octpath_cover_base = octpath_cover.copy()
    # # octpath_cover_base[:,depthBase::]=1
    # # gridlist_cover_base = improc.oct2grid(octpath_cover_base)
    #
    #             # crops
    #             mylist[idTile]
    #
            # TT = dset[0:600, 3712:3712+leafSize[1], 988:988+leafSize[2], 0]
            # cv = np.max(TT, axis=2)
            #
            # fig = plt.figure()
            #     plt.imshow(cv.T)
    #
    #             # xyztest = xres[np.where(all(octpath == octpath[0, :], axis=1))[0], :]
    #             # plt.scatter(xyztest[:, 0], xyztest[:, 1], c='r')
    #
    #
    #
    #
    #
    #
    #             # '/'.join(a + b for a, b in zip(s[::1], s[1::1]))
    #
    #
    # # -> search upto depth to find unique tiles
    # baselist = octpath_dilated[:, :depthBase.__int__()]
    # basetiles = np.unique(baselist, axis=0)
    #
    #
    #
    #
    #
    # listTiles(octpath_dilated,params["nlevels"].__int__())
    #
    #
    #
    #
    # # find bbox for neuron, find envelop
    # BBox = findBBox(octpath_dilated)
    #
    #
    #
    # # TODO:
    # # -> appen to list
    #
    #
    #
    # # fd = dist.cdist(octpath_cover, octpath_dilated)
    # # find bounding box that is rounded to octree format
    # # dilate octree with 1
    # octpath_dilated = improc.dilateOct(octpath_cover)
    # # bounding box
    # gridxyz=improc.oct2grid(octpath_dilated)
    # bbox = np.stack((np.min(gridxyz,axis=0),np.max(gridxyz,axis=0)),axis=0)
    #
    # ##
    # #load tif stack
    # im = io.imread('/nrs/mouselight/SAMPLES/2017-06-10/2/1/3/8/6/1/default.0.tif')
    # imp = np.max(im, axis=0)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(xyzup[:,0], xyzup[:,1], xyzup[:,2],marker='o',c='b')
    # ax.set_xlim([11000, 12000])
    # ax.set_ylim([3700, 3800])
