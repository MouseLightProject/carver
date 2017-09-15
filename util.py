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

def um2pix(um,A):
    # applies transform to convert um into pix location
    return(np.dot(np.linalg.pinv(A),um))
def pix2um(xyz,A):
    return(np.dot(A,xyz))




