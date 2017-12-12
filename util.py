import numpy as np
from collections import defaultdict

def readTransfrom(transformfile = "/nrs/mouselight/SAMPLES/2017-06-10/transform.txt"):
    # reads transform.txt file and parse it into a transform
    A = np.zeros((3,4))
    with open(transformfile, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            parts = text.split(':')
            num = np.float(parts[1].strip('\n'))
            if parts[0] == 'ox':
                A[0,3] = num
            elif parts[0] == 'oy':
                A[1, 3] = num
            elif parts[0] == 'oz':
                A[2,3] = num
            elif parts[0] == 'sx':
                A[0, 0] = num
            elif parts[0] == 'sy':
                A[1, 1] = num
            elif parts[0] == 'sz':
                A[2, 2] = num
            elif parts[0] == 'nl':
                # normalize diagonal with level
                np.fill_diagonal(A,A.diagonal()/(2**(num-1)))
    return A

def readParameterFile(parameterfile = "/nrs/mouselight/SAMPLES/2017-06-10/calculated_parameters.jl"):
    # reads calculated_parameters.txt file and parse it into a transform
    # const jobname = "ocJHfFH"
    # const nlevels = 6
    # const nchannels = 2
    # const shape_leaf_px = [406,256,152]
    # const voxelsize_used_um = [0.329714,0.342888,1.00128]
    # const origin_nm = [61677816,45421726,16585827]
    # const tile_type = convert(Cint,1)
    # const render_version = "2017-03-06 14:11:38 -0500 ec13bbfa7f9285447d3b9702b96a1f1afb847244"
    # const mltk_bary_version = "2016-11-11 12:16:03 -0500 84e153640047e3830abf835e1da4b738efa679d3"
    # const tilebase_version = "2016-08-22 15:49:39 -0400 cc171869a904e9e876426f2bb2732a38e607a102"
    # const nd_version = "2016-11-17 08:30:01 -0500 ef4923831c7bddadd0bba6b706f562a7cde00183"
    # const ndio_series_version = "2016-08-23 11:11:13 -0400 fdfe30a71f3d97fad6ac9982be50d8aea90b5234"
    # const ndio_tiff_version = "2016-08-23 11:11:54 -0400 df46d485cdf55ba66b8ed16fcf9fd9f3d5892464"
    # const ndio_hdf5_version = "2016-08-30 14:25:54 -0400 0c7ac77c5ca535913bfae5300159e6bdf60e36ca"
    # const mylib_version = "2013-08-06 19:15:35 -0400 0ca27aae55a5bab44263ad2e310e8f4507593ddc"
    params = {} # initialize dictionary
    with open(parameterfile, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            parts = text.split('=')
            keyval = parts[0].strip()
            if keyval == 'const nlevels':
                params['nlevels'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const shape_leaf_px':
                params['leafshape'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const voxelsize_used_um':
                params['vixsize'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const origin_nm':
                params['origin'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            elif keyval == 'const nchannels':
                params['nchannels'] = np.array(eval(parts[1].strip('\n')),dtype=np.float)
            else:
                it=0
    A = np.zeros((3,4))
    np.fill_diagonal(A, params['vixsize']*1000) #convert to nm
    A[:,3] = params['origin']
    params['A'] = A
    return params

def readSWC(swcfile='/groups/mousebrainmicro/mousebrainmicro/cluster/Tracings/2017-06-10_G-029_Consensus.swc',scale=1):
    swcline=[]
    offset = np.zeros((1,3))
    offsetkey = 'OFFSET'
    with open(swcfile, 'r') as f:
        while True:
            text = f.readline()
            if not text: break
            if text[0]=='#':
                # check offset
                if text[2:len(offsetkey)+2]==offsetkey:
                    offset = np.array(text[len(offsetkey) + 3:-1].split(), dtype=np.float).reshape(1,3)
                else:
                    continue #header
            else:
                parts = text.split(' ')
                swcline.append(parts)
    lines = np.array(swcline, dtype=float).reshape(-1, 7)
    edges = lines[:,(0,6)]
    R = lines[:,5]
    xyz = lines[:,2:5]
    xyz = xyz + offset
    xyz = xyz/scale
    return (xyz,edges,R,offset,scale)

def upsampleSWC(xyz,edges,sp):
    xyzup = []
    for i, j in np.asarray(edges - 1, np.int64):
        if j < 0:
            continue
        else:
            st = xyz[i, :][None, :]
            ed = xyz[j, :][None, :]
            el = ed - st
            enel = el / np.linalg.norm(el)
            numiter = np.ceil(np.linalg.norm(el) / sp)
            xyzup.append(np.arange(0, numiter).reshape(-1, 1) * enel * sp + st)

    return(np.concatenate(xyzup))


def um2pix(um,A):
    # applies transform to convert um into pix location
    # um_ = np.concatenate((um, np.ones((um.shape[0], 1))), axis=1)
    # return(np.dot(np.linalg.pinv(A),um.T))
    return(np.dot(np.diag(1 / np.diagonal(A[:3, :3])), (um - A[:, 3]).T))

def pix2um(xyz,A):
    return(np.dot(A,xyz))
def pix2oct(xyz,dims,depth):
    # for a given xyz, box size and depth, returns the location int the patch and patch path
    res = dims/depth
    ijk = np.floor(xyz/res)
    # convert ijk to

    return 0

def um2oct(xyz,dims,transform ):
    # for a given um, transform and image size, returns the patch location
    return 0
def traverseOct():
    # lets you to traverse octree
    return 0




