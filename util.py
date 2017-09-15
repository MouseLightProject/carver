import numpy as np
def readTransfrom(transformfile = "/nrs/mouselight/SAMPLES/2017-06-10/transform.txt"):
    # reads transform.txt file and parse it into a transform
    A = np.zeros((3,4))
    A.diagonal()

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

def um2pix(um,A):
    # applies transform to convert um into pix location
    # um_ = np.concatenate((um, np.ones((um.shape[0], 1))), axis=1)
    # return(np.dot(np.linalg.pinv(A),um.T))
    return(np.dot(np.diag(1 / np.diagonal(A[:3, :3])), (um - A[:, 3]).T))
def pix2um(xyz,A):
    return(np.dot(A,xyz))




