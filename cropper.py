import improc
import os
import util
import re
import numpy as np
import h5py
import skimage.io as io


def crop_from_render(data_fold,input_swc,output_folder,output_swc_name,output_h5_name,cast2vox=False):
    output_swc_file = os.path.join(output_folder,output_swc_name)
    output_h5_file =  os.path.join(output_folder,output_h5_name)

    params = util.readParameterFile(parameterfile=data_fold+"/calculated_parameters.jl")
    scale = 1 # for voxel, use 1, for scope use 1/1000 to cast to nm
    nm, edges, R, offset, scale, header = util.readSWC(swcfile=input_swc,scale=scale)

    if cast2vox:
        # to fix the bug in Janelia Workstation
        nm = nm + params['vixsize']/scale/2
        xyz = util.um2pix(nm,params['A']).T
        # LVV BUG FIX: if source is JW, fix coordinates xyz_correct = xyz_LVV-[1 1 0]
        # if any([re.findall('Janelia Workstation Large Volume Viewer', lines) for lines in header]):
        #     xyz=xyz-[1,1,0]

        # upsample xyz to
        xyz = util.upsampleSWC(xyz, edges, sp=5)
    else:
        xyz = nm

    if False:
        octpath, xres = improc.xyz2oct(xyz,params)
    else:
        depthextend = 3
        params_p1=params.copy()
        params_p1["nlevels"]=params_p1["nlevels"]+depthextend
        params_p1["leafshape"]=params_p1["leafshape"]/(2**depthextend)
        octpath, xres = improc.xyz2oct(xyz,params_p1)

    depthBase = params["nlevels"].astype(int)
    depthFull = params_p1["nlevels"].astype(int)
    tileSize = params["leafshape"].astype(int)
    leafSize = params_p1["leafshape"].astype(int)

    octpath_cover = np.unique(octpath, axis=0)
    gridlist_cover = improc.oct2grid(octpath_cover)
    octpath_dilated,gridlist_dilated = improc.dilateOct(octpath_cover)
    #second pass (somewhat heuristic, helps with cropping later on)
    # octpath_dilated,gridlist_dilated = improc.dilateOct(octpath_dilated)

    tilelist = improc.chunklist(octpath_dilated,depthBase) #1..8
    tileids = list(tilelist.keys())

    gridReference = np.min(gridlist_dilated, axis=0)
    volReference = gridReference*leafSize
    gridSize = np.max(gridlist_dilated, axis=0) - gridReference +1
    outVolumeSize = np.append(gridSize*leafSize,2) #append color channel
    chunksize = np.append(leafSize,2)

    # save into cropped swc
    xyz_shifted = xyz-volReference
    with open(output_swc_file,'w') as fswc:
        for iter,txt in enumerate(xyz_shifted):
            fswc.write('{:.0f} {:.0f} {:.4f} {:.4f} {:.4f} {:.2f} {:.0f}\n'.format(edges[iter,0].__int__(),1,txt[0]-1,txt[1],txt[2],1,edges[iter,1].__int__()))

    # write into h5
    with h5py.File(output_h5_file, "w") as f:
        dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0],7), dtype='f')
        for iter, xyz_ in enumerate(xyz_shifted):
            dset_swc[iter,:] = np.array([edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])

        dset = f.create_dataset("volume", outVolumeSize, dtype='uint16', chunks=tuple(chunksize), compression="gzip", compression_opts=9)
        # crop chuncks from a tile read in tilelist
        for iter,idTile in enumerate(tileids):
            print('{} : {} out of {}'.format(idTile, iter, len(tileids)))
            tilename = '/'.join(a for a in idTile)
            tilepath = data_fold+'/'+tilename

            ijkTile = np.array(list(idTile), dtype=int)
            xyzTile = improc.oct2grid(ijkTile.reshape(1, len(ijkTile)))
            locTile = xyzTile * tileSize
            locShift = np.asarray(locTile - volReference,dtype=int).flatten()
            if os.path.isdir(tilepath):

                im = improc.loadTiles(tilepath)
                relativeDepth = depthFull - depthBase

                # patches in idTiled
                for patch in tilelist[idTile]:
                    ijk = np.array(list(patch),dtype=int)
                    xyz = improc.oct2grid(ijk.reshape(1, len(ijk))) # in 0 base

                    start = np.ndarray.flatten(xyz*leafSize)
                    end = np.ndarray.flatten(start + leafSize)
                    # print(start,end)
                    imBatch = im[start[0]:end[0],start[1]:end[1],start[2]:end[2],:]

                    start = start + locShift
                    end = end + locShift
                    dset[start[0]:end[0],start[1]:end[1],start[2]:end[2],:] = imBatch

                    # imgplot = plt.imshow(np.max(im[..., 0], axis=2))

                    # # convert to tif
    # with h5py.File(output_h5_file, "r") as f:
    #     dset = f['volume']
    #     io.imsave(cropped_tif_file, np.swapaxes(dset[:,:,:,0],2,0))


