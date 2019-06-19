import improc
import os
import util
import re
import numpy as np
import h5py
import skimage.io as io

def crop_from_render(data_fold,input_swc,output_folder,output_swc_name,output_h5_name,scale,cast2vox=True):
    output_swc_file = os.path.join(output_folder,output_swc_name)
    output_h5_file =  os.path.join(output_folder,output_h5_name)

    params = util.readParameterFile(parameterfile=data_fold+"/calculated_parameters.jl")
    depthBase = params["nlevels"].astype(int)
    tileSize = params["leafshape"].astype(int)

    # check if input argument is file or folder
    if os.path.isfile(input_swc):
        inputfolder, swc_name_w_ext = os.path.split(input_swc)
        nm, edges, R, offset, scale, header = util.readSWC(swcfile=os.path.join(inputfolder, swc_name_w_ext),
                                                           scale=scale)
    elif os.path.isdir(input_swc):
        inputfolder = input_swc
        nm, edges, R = util.appendSWCfolder(inputfolder,scale=scale) # somewhat redundant but cleaner
        nm_, edges_, R_, filenames, header = util.readSWCfolder(inputfolder,scale=scale)

    if cast2vox:
        # to fix the bug in Janelia Workstation
        nm = nm + params['vixsize']/scale/2
        xyz_ori = util.um2pix(nm,params['A']).T
        # LVV BUG FIX: if source is JW, fix coordinates xyz_correct = xyz_LVV-[1 1 0]
        # if any([re.findall('Janelia Workstation Large Volume Viewer', lines) for lines in header]):
        #     xyz=xyz-[1,1,0]
        # upsample xyz if needed
        centerpoint_count = xyz_ori.shape[0]
        if centerpoint_count > 1 :
            diff_xyz = np.diff(xyz_ori,axis=0)
            norm_xyz = np.sqrt(np.sum(diff_xyz**2,axis=1))
            if np.mean(norm_xyz)>5:
                xyz = util.upsampleSWC(xyz_ori, edges, sp=10)
            else:
                xyz = xyz_ori
        else:
            xyz = xyz_ori
    else:
        xyz = nm

    # Each tile in the rendered image will itself be 'octreed' into a set of 'leafs'
    # depthextend tells now many octree levels there will be within each tile

    depthextend = 3
    params_p1 = params.copy()
    params_p1["nlevels"]=params_p1["nlevels"]+depthextend
    params_p1["leafshape"]=params_p1["leafshape"]/(2**depthextend)
    octpath, xres = improc.xyz2oct(xyz,params_p1)

    depthFull = params_p1["nlevels"].astype(int)
    leafSize = params_p1["leafshape"].astype(int)

    octpath_cover = np.unique(octpath, axis=0)
    gridlist_cover = improc.oct2grid(octpath_cover)

    octpath_dilated = octpath_cover
    desired_carve_out_half_diagonal_as_scalar = 512
    desired_carve_out_half_diagonal = desired_carve_out_half_diagonal_as_scalar * np.array([1.0, 1.0, 1.0/4.0])
    #dilation_count = 8
    dilation_count = np.max( np.ceil(desired_carve_out_half_diagonal.astype(float) / leafSize.astype(float)) ).astype(int)
    # should be enough to get about a 512 vx cube around each swc centerpoint
    # (except 4x less in z, b/c axial rez is less)
    for dilation_index in range(dilation_count):
        octpath_dilated, gridlist_dilated = improc.dilateOct(octpath_dilated)


    tilelist = improc.chunklist(octpath_dilated,depthBase) #1..8

    tileids = list(tilelist.keys())
    # base on bounding box (results in cropped output volume)
    # gridReference = np.min(gridlist_dilated, axis=0)
    # gridSize = np.max(gridlist_dilated, axis=0) - gridReference +1
    # base on initial image
    gridReference = np.array((0,0,0))
    gridSize = tileSize*(2**(depthBase))/leafSize

    volReference = gridReference*leafSize
    outVolumeSize = np.append(gridSize*leafSize,2) #append color channel
    chunksize = np.append(leafSize,2)

    setting = dict()
    setting['volSize'] = outVolumeSize
    setting['chunkSize'] = tuple(chunksize)
    setting['depthBase'] = depthBase
    setting['depthFull'] = depthFull
    setting['tileSize'] = tileSize
    setting['leafSize'] = leafSize
    setting['volReference'] = volReference

    setting['compression'] = "gzip"
    setting['compression_opts'] = 9
    setting['dtype'] = 'uint16'
    setting['type'] = 'h5'

    # save into cropped swc
    xyz_shifted = xyz_ori-volReference
    if xyz_shifted.ndim<2:
        xyz_shifted = xyz_shifted[None,:] # extend dimension to prevent r < h cases

    # TODO: fix edge upsampling before dumping swc file
    # with open(output_swc_file,'w') as fswc:
    #     for iter,txt in enumerate(xyz_shidefaultDataPlaceholderfted):
    #         fswc.write('{:.0f} {:.0f} {:.4f} {:.4f} {:.4f} {:.2f} {:.0f}\n'.format(edges[iter,0].__int__(),1,txt[0]-1,txt[1],txt[2],1,edges[iter,1].__int__()))

    # dump into file
    with h5py.File(output_h5_file, "w") as f:
        dset_swc = f.create_dataset("reconstruction", (xyz_shifted.shape[0], 7), dtype='f')
        for iter, xyz_ in enumerate(xyz_shifted):
            dset_swc[iter, :] = np.array(
                [edges[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges[iter, 1].__int__()])

    if os.path.isdir(input_swc):
        # dump into file
        with h5py.File(output_h5_file, "a") as f:
            swc_group = f.create_group("swc_files")

            f_swc_original = swc_group.create_group("original")
            for it, swcname in enumerate(filenames):
                numrows = R_[it].shape[0]
                dset_swc = f_swc_original.create_dataset(swcname,(numrows, 7), dtype='f')
                edges_swc = edges_[it]
                nm_swc = nm_[it]
                for iter, xyz_ in enumerate(nm_swc):
                    dset_swc[iter, :] = np.array(
                        [edges_swc[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges_swc[iter, 1].__int__()])

            f_swc_shifted = swc_group.create_group("shifted")
            for it, swcname in enumerate(filenames):
                numrows = R_[it].shape[0]
                dset_swc = f_swc_shifted.create_dataset(swcname,(numrows, 7), dtype='f')
                edges_swc = edges_[it]
                nm_swc = nm_[it]
                nm_swc = nm_swc + params['vixsize'] / scale / 2
                xyz_ori = util.um2pix(nm_swc, params['A']).T
                nm_swc = xyz_ori - volReference

                for iter, xyz_ in enumerate(nm_swc):
                    dset_swc[iter, :] = np.array(
                        [edges_swc[iter, 0].__int__(), 1, xyz_[0], xyz_[1], xyz_[2], 1.0, edges_swc[iter, 1].__int__()])

    # Finally, write the carved data to disk
    #dump = util.dumper(data_fold, output_h5_file, setting,tilelist=tilelist)
    #dump.write()
    util.dump_write(data_fold, output_h5_file, setting, tilelist)

# end def crop_from_render()
