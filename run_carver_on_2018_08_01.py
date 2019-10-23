#! /usr/bin/env python

import carver

sample_date = '2018-08-01'

output_volume_path = '/nrs/funke/mouselight-v3/%s/fluorescence-near-consensus.n5' % sample_date
render_folder_path = '/nrs/mouselight/SAMPLES/%s' % sample_date
swc_folder_path = '/groups/mousebrainmicro/mousebrainmicro/tracing_complete/%s' % sample_date

# First pass is parallel, does most of the work, but usually errors out for some reason
do_use_simple_for_loop = False
try:
    carver.crop_from_render(output_volume_path,
                            render_folder_path,
                            swc_folder_path,
                            do_use_simple_for_loop)
except RuntimeError:
    print("Some dumb Dask runtime error happened")

# Do a serial pass to finish up
do_use_simple_for_loop = True
carver.crop_from_render(output_volume_path,
                        render_folder_path,
                        swc_folder_path,
                        do_use_simple_for_loop)
