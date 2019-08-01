import numpy as np
import math
import re
from collections import defaultdict
import os
from skimage import io
import tqdm



def boundingbox(xyz):
    # finds the bounding box of point cloud
    return(np.round(np.array([np.min(xyz,axis=1),np.max(xyz,axis=1)])))

def snapLoc(I,nploc,w=1):
    # finds the location of maxima in a 3x3x3 neighborhood
    nplocupdated = nploc.copy()
    for iter,loc in enumerate(nploc):
        Icrop = I[loc[0]-w:loc[0]+w+1,
                  loc[1]-w:loc[1]+w+1,
                  loc[2]-w:loc[2]+w+1]
        i, j, k = np.unravel_index(Icrop.argmax(), Icrop.shape)
        nplocupdated[iter] = loc+(np.array((i,j,k))-w)
    return nplocupdated


def ijk2oct(ijk, level_count, chunk_shape):
    # converts ijk location to oct location
    # ALT: params['nlevels'] is the number of levels in the octree, i.e. the length of the path to a leaf stack
    # ALT: params['leafshape'] is the shape of the leaf stacks in the octree, in ijk/xyz order
    # ALT: ijk is in voxels, and is zero-based.  Values should be integral, although the type need not be, I don't think.
    # ALT" ijk is n x 3
    # ALT: chunk_shape is the chunk shape at the bottom level of the tree, i.e. the most zoomed-in level
    # ALT: Returns two things: First is the octree path, an n x nlevels array, row i the path to the leaf stack for ijk[i]
    # ALT:                     Second n x 3, row i giving the ijk position in the leaf stack for ijk[i]
    if ijk.ndim==1:
        did_promote = True
        ijk = ijk[None,:]  # convert vector to a single-row matrix
    else:
        did_promote = False

    level_count = np.int(level_count)
    #chunk_shape = params['leafshape']
    row_count = ijk.shape[0]
    octpath = np.zeros((row_count,level_count), dtype='int64')
    ijk_within_chunk = np.zeros((row_count,3), dtype='int64')
    for idx in range(row_count):
        bits = []
        ijk_this = ijk[idx]
        #u = chunk_shape
        for n in range(level_count-1,-1,-1):
            bn = (2**n)*chunk_shape
            th = ijk_this>bn
            bits.append(th)
            ijk_this = ijk_this - bn*th
        # convert to octodigit
        octpath[idx,:] = (1 + np.sum(np.array(bits)*2**np.array([0,1,2]),axis=1))[None,:]
        ijk_within_chunk[idx,:] = ijk_this

    if did_promote:
        # demote back to vectors
        octpath = np.ndarray.flatten(octpath)
        ijk_within_chunk = np.ndarray.flatten(ijk_within_chunk)

    return octpath, ijk_within_chunk


# def xyz2oct(xyz_in_voxels, level_count, leaf_shape):
#     # converts xyz location to oct location
#     # ALT: params['nlevels'] is the number of levels in the octree, i.e. the length of the path to a leaf stack
#     # ALT: params['leafshape'] is the shape of the leaf stacks in the octree, in xyz order
#     # ALT: xyz is in voxels, and is zero-based
#     # ALT" xyz is n x 3
#     # ALT: Returns two things: First is the octree path, an n x nlevels array, row i the path to the leaf stack for xyz[i]
#     # ALT:                     Second n x 3, row i giving the xyz position in the leaf stack for xyz[i]
#
#     # If xyz_in_voxels is a vector, convert to a one-row matrix
#     if len(xyz_in_voxels.shape)==1:
#         xyz_in_voxels= xyz_in_voxels[None, :]
#
#     level_count = np.int(level_count)
#     row_count = xyz_in_voxels.shape[0]
#     octpath = np.zeros((row_count, level_count))
#     xres = np.zeros((row_count, 3))
#     for row_index in range(row_count):
#         bits = []
#         x = xyz_in_voxels[row_index]
#         u = leaf_shape
#         for n in range(level_count - 1, -1, -1):
#             bn = (2**n)*u
#             th = x>bn
#             bits.append(th)
#             x = x - bn*th
#         # convert to octodigit
#         octpath[row_index,:] = (1+ np.sum(np.array(bits)*2**np.array([0,1,2]),axis=1))[None,:]
#         xres[row_index,:] = x
#
#     return octpath.astype(int),xres


def to_base_3(n):
    s = ""
    while n:
        print(n)
        s = str(n % 3) + s
        n = round(n/3)
    return s

def to_base_2(n,numdigit=0):
    n = math.floor(n)
    s = ""
    while n:
        s = str(n % 2) + s
        n = math.floor(n/2)
    s=(numdigit-len(s))*'0'+s
    return s

def oct2grid_list(octpath):
    depth = octpath.shape[1]
    numpath = octpath.shape[0]
    alltiles = []
    for ijk in octpath:
        xyz = oct2grid(ijk.reshape(1,depth))
        # for every path, there are 26 neighbors
        alltiles.append(xyz[None,:])
    if len(alltiles) == 1:
        alltiles = np.squeeze(alltiles[0])
    else:
        alltiles = np.squeeze(np.concatenate(alltiles, axis=1))
    return alltiles

def oct2grid(oct_idx):
    # (inverse logic as grid2oct)
    # oct_idx [1..8]
    # grid [0 dims]
    if np.any(oct_idx < 1) or np.any(oct_idx > 8):
        raise Exception('oct out of bound, oct \in [1...8]')

    if oct_idx.ndim == 1:
        oct_idx = oct_idx.reshape(1,len(oct_idx))

    numlist = oct_idx.shape[0]
    depth = oct_idx.shape[1]
    binarray = 2 ** (np.array(range(depth, 0, -1)) - 1)
    gridarray = np.zeros((numlist,3))
    for il in range(numlist):
        idxarray = np.zeros((3, depth))
        for id in range(depth):
            base2 = to_base_2(oct_idx[il,id]-1, 3)
            idxarray[2, id] = int(base2[0])
            idxarray[1, id] = int(base2[1])
            idxarray[0, id] = int(base2[2])
        gridarray[il,:] = np.sum(idxarray * binarray, axis=1)
    # broadcast binarray
    return(np.asarray(gridarray,dtype=int))

def loadTiles(tilepath,ext=".tif"):
    IM=[]
    files = os.listdir(tilepath)
    files.sort() # make sure that channels are loaded in order
    for file in files:
        if file.endswith(ext):
            tilefiles = os.path.join(tilepath, file)
            # load tile if exists
            if os.path.isfile(tilefiles):
                im = io.imread(tilefiles)  # zyx order
                IM.append(np.swapaxes(im, 0, 2))

    return np.stack(IM,axis=3)



def grid2oct(ijks, level_count):
    # Converts an n x 3 array of chunk coordinates (in xyz order, each row a set of coordinates) to
    # an n x level_count list of octree paths, as used for mouselight octrees.
    chunk_count = ijks.shape[0]
    zero_based_chunk_paths = np.zeros((chunk_count, level_count), dtype=np.int)
    for il in range(chunk_count):
        arr = []
        arr.append(to_base_2(ijks[il, 2], level_count))
        arr.append(to_base_2(ijks[il, 1], level_count))
        arr.append(to_base_2(ijks[il, 0], level_count))
        for idx in range(level_count):
            b = [el[idx] for el in arr]
            zero_based_chunk_paths[il,idx]=np.int(''.join(b),2)
    chunk_paths = zero_based_chunk_paths+1
    return chunk_paths



def chunklist(pathlist,depth):
    # -> for each tile, find bbox of crop sub-octtree
    listdict = defaultdict(list)
    for tileid in pathlist:
        # list of crops for this tile
        mykey = re.sub('[\[\]]', '', np.array_str(tileid).replace(' ',''))
        listdict[mykey[:depth]].append(mykey[depth:])

    return listdict



def dilateOct(octpath, width=1):
    # dilates the octpath with the given search widty
    # 1/2/3 with width 1 -> 1/2/3 | 2/2/3 | 1/1/3 | 1/3/3 | ... | 2/3/4
    depth = octpath.shape[1]
    # numpath = octpath.shape[0]
    ix, iy, iz = np.mgrid[-width:width + 1, -width:width + 1, -width:width + 1]
    ixyz = np.stack((ix.flatten(), iy.flatten(), iz.flatten()), axis=1)
    alltiles = []
    for ijk in octpath:
        xyz = oct2grid(ijk.reshape(1, depth))
        # for every path, there are 26 neighbors
        alltiles.append(xyz[None, :] + ixyz)

    if len(alltiles) == 1:
        alltiles = np.squeeze(alltiles[0])
    else:
        alltiles = np.squeeze(np.concatenate(alltiles, axis=1))

    # delete any out of bound tiles
    deletethese = np.any(np.logical_or(alltiles < 0, alltiles > 2 ** depth - 1), axis=1)
    alltiles = np.delete(alltiles, (np.where(deletethese)), axis=0)

    # unique entries
    alltiles_unique = np.unique(alltiles, axis=0)

    # convert to octpaths
    octlist = [grid2oct(tileid[None, :], depth) for tileid in alltiles_unique]
    octlist = np.concatenate(octlist, axis=0)

    return octlist



def dilate_octree_chunk_set(octree_chunk_paths, radius):
    # Given a list of paths to octree chunks, determines all the octree chunks within
    # radius using city-block distance, and returns a list of these chunks.

    # Get dimensions
    path_count = octree_chunk_paths.shape[0]
    level_count = octree_chunk_paths.shape[1]
    chunk_stack_size = 2 ** level_count  # the stack of chunks is this x this x this

    # Generate a list of offsets that constitute the neighborhood of a chunk
    delta_i_grid, delta_j_grid, delta_k_grid = np.mgrid[-radius:radius + 1, -radius:radius + 1, -radius:radius + 1]
    neighborhood_delta_ijks = np.stack((delta_i_grid.flatten(), delta_j_grid.flatten(), delta_k_grid.flatten()), axis=1)

    # Mark all the chunks that will be included in the output list
    is_chunk_marked = np.zeros((chunk_stack_size, chunk_stack_size, chunk_stack_size), dtype=bool)
    for octree_chunk_path in tqdm.tqdm(octree_chunk_paths):
        central_ijk = oct2grid(octree_chunk_path.reshape(1,level_count))
        neighborhood_ijks = central_ijk + neighborhood_delta_ijks
        #for neighborhood_ijk in neighborhood_ijks:
        #    i = neighborhood_ijk[0]
        #    j = neighborhood_ijk[1]
        #    k = neighborhood_ijk[2]
        #    if (0<=i) and (i<chunk_stack_size) and (0<=j) and (j<chunk_stack_size) and (0<=k) and (k<chunk_stack_size):
        #        is_chunk_marked[i,j,k] = True
        is_neighborhood_chunk_in_bounds = np.all(
            np.logical_and(0 <= neighborhood_ijks, neighborhood_ijks < chunk_stack_size), axis=1)
        in_bounds_neighborhood_ijks = neighborhood_ijks[is_neighborhood_chunk_in_bounds]
        is_chunk_marked[tuple(np.transpose(in_bounds_neighborhood_ijks))] = True

    # Find all the marked chunks, convert to chunk paths, and return
    ijks_in_output = np.transpose(np.nonzero(is_chunk_marked))  # marked_chunk_count x 3
    result = grid2oct(ijks_in_output, level_count)
    return result
# end


