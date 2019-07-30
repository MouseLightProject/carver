#! /usr/bin/env python

import navigator

navigator.main(['-i', '/nrs/mouselight/SAMPLES/2018-10-01',
                '-s', '/groups/mousebrainmicro/mousebrainmicro/scripts/navigator/2018-10-01/consensus-neurons-with-machine-centerpoints-labelled-as-swcs',
                '-o', '/nrs/mouselight/cluster/navigator-output/2018-10-01-tile-chunks-on-cluster'])
