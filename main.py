import util
import numpy as np
import matplotlib.pyplot as plt
import improc
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

    epsball = 500

    import importlib

    importlib.reload(util)
    importlib.reload(improc)

    params = util.readParameterFile(parameterfile="./calculated_parameters.jl")
    nm, edges, R, offset, scale = util.readSWC(scale=1.0/1000)
    xyz = util.um2pix(nm,params['A']).T
    # # upsample xyz to
    # sp = 5
    # xyzup = util.upsampleSWC(xyz, edges, sp)
    if False:
        octpath, xres = improc.xyz2oct(xyz,params)
    else:
        params_p1=params.copy()
        params_p1["nlevels"]=params_p1["nlevels"]+1
        params_p1["leafshape"]=params_p1["leafshape"]/2
        octpath, xres = improc.xyz2oct(xyz,params_p1)

    octpath_cover = np.unique(octpath, axis=0)
    octpath_dilated = improc.dilateOct(octpath_cover)

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
