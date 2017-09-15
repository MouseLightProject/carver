import util

def crop_from_render():
    inputtif = "/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/compression/compression/testin.tif"
    mj2file = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/compression/compression/output5.avi'
    inputmj = '/groups/mousebrainmicro/home/base/CODE/forNelsonClassfier/outdata/compressed_mj2/granule-1_comp-10.mj2'
    decompression(inputmj)
    # decomp
    # tif2mj2(inputtif,mj2file)

if __name__ == '__main__':
    crop_from_render()