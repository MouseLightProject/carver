import util
import numpy as np
import matplotlib.pyplot as plt
import improc
from scipy.spatial import distance as dist

from mpl_toolkits.mplot3d import Axes3D
from skimage import io

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def crop_from_render():
    # decomp
    # tif2mj2(inputtif,mj2file)
    test = 1

if __name__ == '__main__':

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
        params_p1=params.copy()
        params_p1["nlevels"]=params_p1["nlevels"]+1
        params_p1["leafshape"]=params_p1["leafshape"]/2
        octpath, xres = improc.xyz2oct(xyzup,params_p1)

    octpath_cover = np.unique(octpath, axis=0)
    octpath_dilated = improc.dilateOct(octpath_cover)
    octpath_dilated = np.unique(octpath_dilated, axis=0)

    # TODO:
    # -> search upto depth to find unique tiles
    # -> for each tile, find bbox of crop sub-octtree
    # -> appen to list



    # fd = dist.cdist(octpath_cover, octpath_dilated)
    # find bounding box that is rounded to octree format
    # dilate octree with 1
    octpath_dilated = improc.dilateOct(octpath_cover)
    # bounding box
    gridxyz=improc.oct2grid(octpath_dilated)
    bbox = np.stack((np.min(gridxyz,axis=0),np.max(gridxyz,axis=0)),axis=0)

    xyztest = xres[np.where(np.all(octpath == octpath[0, :], axis=1))[0], :]
    aha = np.asarray(np.round(xyztest[1])[::-1],np.int)

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
