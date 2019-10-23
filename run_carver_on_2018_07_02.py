#! /usr/bin/env python

import carver

sample_date = '2018-07-02'

output_volume_path = '/nrs/funke/mouselight-v3/%s/fluorescence-near-consensus.n5' % sample_date
render_folder_path = '/nrs/mouselight/SAMPLES/%s' % sample_date
swc_folder_path = '/groups/mousebrainmicro/mousebrainmicro/tracing_complete/%s' % sample_date
do_use_simple_for_loop = False

carver.crop_from_render(output_volume_path,
                        render_folder_path,
                        swc_folder_path,
                        do_use_simple_for_loop)

# Do a serial pass to finish up
do_use_simple_for_loop = True
carver.crop_from_render(output_volume_path,
                        render_folder_path,
                        swc_folder_path,
                        do_use_simple_for_loop)
