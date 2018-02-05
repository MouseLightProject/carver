import util
import numpy as np
import matplotlib.pyplot as plt
import improc
import cropper

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

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


