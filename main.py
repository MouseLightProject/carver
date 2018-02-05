import util
import numpy as np
import matplotlib.pyplot as plt
import improc
import cropper


from scipy.spatial import distance as dist
from morphsnakes import tests
from morphsnakes import morphsnakes as morph


from mpl_toolkits.mplot3d import Axes3D

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def snake(img):
    tests.test_confocal3d()

    macwe = morph.MorphACWE(img, smoothing=1, lambda1=1, lambda2=2)
    macwe.levelset = morph.circle_levelset(img.shape, (30, 50, 80), 25)

    viz = 1;
    if viz:
        fig = plt.figure(frameon=False)
    for i in range(100):
        macwe.step()
        if viz:
            fd = macwe.levelset
            plt.imshow(np.max(img, axis=2), cmap='gray')
            plt.imshow(np.max(fd, axis=2), cmap='jet', alpha=0.3)
            plt.show();plt.pause(.1);plt.draw()

if __name__ == '__main__':

    # oct in [1...8]
    # grid in [0...(2**depth-1)]

    import importlib
    importlib.reload(util)
    importlib.reload(improc)

    datafold = '/nrs/mouselight/SAMPLES/2017-11-17'
    outfolder = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data'
    swcfolder = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/recon_repo'
    swcname = '2017-11-17_G-017_Seg-1'
    cropper.crop_from_render(swcfolder,swcname,datafold,outfolder)
    # volumetric segmentation around reconstruction

    tests.test_confocal3d()


