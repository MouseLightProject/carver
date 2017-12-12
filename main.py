import util
import numpy as np
import matplotlib.pyplot as plt
import improc
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
    octpath, xres = improc.xyz2oct(xyz,params)

    # load tile
    kk=0
    sample = '/nrs/mouselight/SAMPLES/2017-06-10/'
    tif0 = sample+'/'.join(np.str(int(a)) for a in octpath[kk])+'/default.0.tif'

    im = io.imread(tif0)
    imp = np.max(im, axis=0)

    ik=0
    loc = np.asarray(xres[0:2,:],dtype=np.int)
    locupdated = improc.snapLoc(im,loc)

    frag_x = locupdated[:,0]
    frag_y = locupdated[:,1]
    frag_z = locupdated[:,2]

    # frag_x = loc[:,0]
    # frag_y = loc[:,1]
    # frag_z = loc[:,2]

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharex=ax1)
    ax1.imshow(imp)
    plt.hold(True)
    ax1.plot(frag_x,frag_y,'r')
    ax2.imshow(np.max(im[np.int(frag_z[0]):np.int(frag_z[1])],axis=0))
    plt.hold(True)
    ax2.plot(frag_x,frag_y,'r')





    # find bounding box that is rounded to octree format
    octpath_cover = np.unique(octpath, axis=0)
    # dilate octree with 1
    octpath_dilated = improc.dilateOct(octpath_cover)
    # bounding box
    gridxyz=improc.oct2grid(octpath_dilated)
    bbox = np.stack((np.min(gridxyz,axis=0),np.max(gridxyz,axis=0)),axis=0)

    xyztest = xres[np.where(np.all(octpath == octpath[0, :], axis=1))[0], :]
    aha = np.asarray(np.round(xyztest[1])[::-1],np.int)



    # plt.imshow(im[aha[0], :, :], interpolation='none')
    # plt.scatter(aha[2], aha[1], c='r')
    # fig = plt.figure()
    # for ii in range(10):
    #     aha = np.asarray(np.round(xyztest[ii])[::-1], np.int) #-1 is heuris
    #     plt.subplot(351+ii)
    #     plt.imshow(im[aha[0],:,:], interpolation='none')
    #     plt.scatter(aha[2],aha[1],c='r')
    #
    #     # fig = plt.figure()
    #     # plt.subplot(131)
    #     # plt.imshow(im[aha[0],:,:])
    #     # plt.scatter(aha[2],aha[1],c='r')
    #     # plt.subplot(132)
    #     # plt.imshow(im[aha[0]+1,:,:])
    #     # plt.scatter(aha[2],aha[1],c='r')
    #     # plt.subplot(133)
    #     # plt.imshow(im[aha[0]+2,:,:])
    #     # plt.scatter(aha[2],aha[1],c='r')
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
