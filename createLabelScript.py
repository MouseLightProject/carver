import util
import numpy as np
import cropper
import os
import h5py
import sys, getopt

input_folder = '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo'
swcfiles = [os.path.join(input_folder, fold, files) for fold in os.listdir(input_folder) if
           os.path.isdir(os.path.join(input_folder, fold)) for files in os.listdir(os.path.join(input_folder, fold)) if
           files.endswith("-carved.swc")]
swcfiles.sort()
for swc_file in swcfiles:
    path, filename = os.path.split(swc_file)
    output_h5_file = os.path.join(path, filename.split('.')[0] + '.h5')
    input_swc = swc_file
    # output_h5_file = os.path.join(path,filename.split('.')[0][:-1]+'.h5')
    # input_swc = os.path.join(path,filename.split('.')[0][:-1]+'.swc')
    with h5py.File(output_h5_file, "r+") as f:
        try:
            del f['reconstruction']
        except Exception:
            pass
        um, edges, R, offset, scale, header = util.readSWC(swcfile=input_swc, scale=1)
        dset_swc = f.create_dataset("reconstruction", (um.shape[0], 7), dtype='f')
        for iter, xyz_ in enumerate(um):
            xyz_ = np.ceil(xyz_-np.sqrt(np.finfo(float).eps))
            dset_swc[iter, :] = np.array([edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])
