import util
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def crop_from_render():
    inputtif = "/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/compression/compression/testin.tif"
    mj2file = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/compression/compression/output5.avi'
    inputmj = '/groups/mousebrainmicro/home/base/CODE/forNelsonClassfier/outdata/compressed_mj2/granule-1_comp-10.mj2'
    decompression(inputmj)
    # decomp
    # tif2mj2(inputtif,mj2file)

if __name__ == '__main__':

    import importlib
    importlib.reload(util)

    A = util.readTransfrom()
    (um, edges, R, offset, scale) = util.readSWC(scale=1/1000)
    # return (np.dot(np.linalg.pinv(A), um.T))
    xyz = util.um2pix(um,A)
    ## plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xyz[0],xyz[1],xyz[2],'r.')





