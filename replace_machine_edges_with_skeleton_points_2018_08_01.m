sample_date = '2018-08-01' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/%s/carver', sample_date) ;
fragments_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/prob0_swcs/frags', sample_date) ;
consensus_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;
sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
skeleton_folder_path =  sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;

replace_machine_edges_with_skeleton_points(output_folder_path, fragments_folder_path, consensus_swcs_folder_path, sample_folder_path, skeleton_folder_path) ;
