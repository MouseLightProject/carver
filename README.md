# navigator

creates cropped volume and JW structure (for visualization) based on input render folder and swc file

# usage: 
```
'navigator.py -i <data_folder> -s <swc_file> -o <output_folder>'  
        -i <data_folder>: input data folder. Folders should follow octree format, e.g. <data_folder>/1/5/6
        -s <swc_file>: input swc_file or folder. for *swc files 7 column conventional reconstruction format.
        -o <output_folder>: folder to create h5 and JW files
        -h <number_of_level>: [OPTIONAL] sets how many chunks around trace will be used
```

## example:
```
python navigator.py -i /nrs/mouselight/SAMPLES/2018-08-01-raw-rerender -s '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data/swc_recons/2018-08-01' -o '/groups/mousebrainmicro/mousebrainmicro/users/base/AnnotationData/h5repo/2018-08-01'
```

# notes:
        Don't use trailing slash '/'. Somehow it can mass with created h5 name. `
        oct in [1...8]
        grid in [0...(2**depth-1)]

        we keep mouselight data in <root>/<neuron-id>/consensus/<tag>_consensus.swc format, e.g.:
        /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-08-01/G-002/consensus/2018-08-01_G-002_consensus.swc
        it is suggested to copy all consensus files for that sample into a single folder manually or with a script than pass input folder with "-f" argument.
        For example:
        cd /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-08-01
        find . -name "*consensus*.swc" -exec cp {} /groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/navigator/data/swc_recons/2018-08-01 \;
