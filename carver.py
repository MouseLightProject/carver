#! /usr/bin/env python

"""
Front-end script for making a carve-out volume
"""

import argparse
import os
import subprocess
import sys
import navigator

parser = argparse.ArgumentParser(description='Carve-out data near consensus neurons for a brain')
parser.add_argument('sample_date', help='Date of the sample to use, e.g. "2019-05-27"')
args = parser.parse_args()

sample_date = args.sample_date
matlab_command_line_template = 'try; label_machine_centerpoints_given_sample_date(\'%s\'); catch err; fprintf(2, \'%%s\\n\', err.getReport()); quit(1); end; quit(0);'
print(matlab_command_line_template)
matlab_command_line = matlab_command_line_template % (sample_date)
print(matlab_command_line)

script_file_path = os.path.abspath(__file__)
script_folder_path = os.path.dirname(script_file_path)
os.chdir(script_folder_path)
child = subprocess.Popen(['/misc/local/matlab-2018b/bin/matlab', '-nodisplay', '-r', matlab_command_line])
child.communicate()
rc = child.returncode
if rc != 0:
    sys.exit(rc)

rendered_sample_path = '/nrs/mouselight/SAMPLES/%s' % (sample_date)
swc_folder_path = '/nrs/mouselight/cluster/classifierOutputs/%s/consensus-neurons-with-machine-centerpoints-labelled/as-swcs' % (sample_date)
output_folder_path = '/nrs/funke/mouselight-v2/%s' % (sample_date)


navigator.main(['-i', rendered_sample_path,
                '-s', swc_folder_path,
                '-o', output_folder_path])
