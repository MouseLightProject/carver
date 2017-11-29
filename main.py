import util
import numpy as np
import matplotlib.pyplot as plt
import improc
import os
from skimage import io
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

    params = util.readParameterFile(parameterfile="./calculated_parameters.jl")
    nm, edges, R, offset, scale = util.readSWC(scale=1.0/1000)
    xyz = util.um2pix(nm,params['A']).T
    # # upsample xyz to
    sp = 10
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
    octpath_cover = np.unique(octpath, axis=0)

    octpath_dilated,gridlist = improc.dilateOct(octpath_cover)
    gridlist_cover = improc.oct2grid_list(octpath_cover)

    # reference shift
    octpath_cover_base = octpath_cover.copy()
    octpath_cover_base[:,depthBase::]=1
    gridlist_cover_base = improc.oct2grid_list(octpath_cover_base)



    gridsize = np.max(gridlist, axis=0) - np.min(gridlist, axis=0) +1
    outvolumeSize = gridsize*params_p1["leafshape"]
    chunksize = params_p1["leafshape"].astype(int)

    # f = h5py.File("test.hdf5", "w")
    # dset = f.create_dataset("volume", outvolumeSize, chunks=tuple(chunksize))
    # f.__delitem__("volume")

    mylist = improc.chunklist(octpath_dilated,depthBase.__int__())
    tileids = list(mylist.keys())


    for file in os.listdir("/mydir"):
        if file.endswith(".txt"):
            print(os.path.join("/mydir", file))


    datafold = '/Volumes/mouselight/SAMPLES/2017-06-10'
    for idTile in tileids:
        tilename = '/'.join(a for a in idTile)
        tilepath = datafold+'/'+tilename
        for file in os.listdir(tilepath):
            if file.endswith(".tif"):
                tilefiles = os.path.join(tilepath, file)
                # load tile
                im = io.imread(tilefiles) # zyx order
                im2=np.swapaxes(im, 0, 2)

                # crops
                mylist[idTile]

                fig = plt.figure()
                plt.imshow(np.max(im2,axis=2))

                # xyztest = xres[np.where(all(octpath == octpath[0, :], axis=1))[0], :]
                # plt.scatter(xyztest[:, 0], xyztest[:, 1], c='r')






                # '/'.join(a + b for a, b in zip(s[::1], s[1::1]))


    # -> search upto depth to find unique tiles
    baselist = octpath_dilated[:, :depthBase.__int__()]
    basetiles = np.unique(baselist, axis=0)





    listTiles(octpath_dilated,params["nlevels"].__int__())




    # find bbox for neuron, find envelop
    BBox = findBBox(octpath_dilated)



    # TODO:
    # -> appen to list



    # fd = dist.cdist(octpath_cover, octpath_dilated)
    # find bounding box that is rounded to octree format
    # dilate octree with 1
    octpath_dilated = improc.dilateOct(octpath_cover)
    # bounding box
    gridxyz=improc.oct2grid(octpath_dilated)
    bbox = np.stack((np.min(gridxyz,axis=0),np.max(gridxyz,axis=0)),axis=0)

    ##
    #load tif stack
    im = io.imread('/nrs/mouselight/SAMPLES/2017-06-10/2/1/3/8/6/1/default.0.tif')
    imp = np.max(im, axis=0)

    fig = plt.figure()
    plt.imshow(imp)
    xyztest = xres[np.where(all(octpath == octpath[0, :], axis=1))[0], :]
    plt.scatter(xyztest[:,0],xyztest[:,1],c='r')
    ##
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],marker='+',c='r')
    # ax.set_xlim([11000, 12000])
    # ax.set_ylim([3700, 3800])
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(xyzup[:,0], xyzup[:,1], xyzup[:,2],marker='o',c='b')
    # ax.set_xlim([11000, 12000])
    # ax.set_ylim([3700, 3800])
