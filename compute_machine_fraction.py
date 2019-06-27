import os
import numpy as np
import util

def read_swc_structure_identifier_column(file_name):
    (xyz,edges,R,offset,scale,header,structure_identifier) = util.readSWC(file_name)
    return structure_identifier

def read_swc_is_centerpoint_machine_generated(file_name):
    structure_identifier_column = read_swc_structure_identifier_column(file_name)
    is_machine_generated = np.equal(structure_identifier_column, 43)
    return is_machine_generated

swc_folder_path = '/groups/scicompsoft/home/taylora/cache/gt/2018-10-01/consensus-neurons-with-machine-centerpoints-labelled-as-swcs'
swc_file_names = os.listdir(swc_folder_path)

swc_file_paths = [os.path.join(swc_folder_path, swc_file_name) for swc_file_name in swc_file_names]
is_machine_generated_from_centerpoint_id_from_file_id = [read_swc_is_centerpoint_machine_generated(swc_file_path) for swc_file_path in swc_file_paths]
machine_generated_centerpoint_count_from_file_id = np.array(list(map(np.sum, is_machine_generated_from_centerpoint_id_from_file_id)))
centerpoint_count_from_file_id = np.array(list(map(np.size, is_machine_generated_from_centerpoint_id_from_file_id)))

fraction_machine_generated_from_file_id = machine_generated_centerpoint_count_from_file_id / centerpoint_count_from_file_id
overall_fraction_machine_generated = np.sum(machine_generated_centerpoint_count_from_file_id) / np.sum(centerpoint_count_from_file_id)
min_fraction_machine_generated = np.min(fraction_machine_generated_from_file_id)
max_fraction_machine_generated = np.max(fraction_machine_generated_from_file_id)

