sample_date = '2018-07-02' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
output_folder_path = fullfile(this_folder_path, sample_date) ;
cache_folder_path = fullfile(this_folder_path, sample_date) ;
fragments_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/prob0_swcs/frags', sample_date) ;
consensus_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;
sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;

consensus_neurons_with_machine_centerpoints_labelled = ...
    label_machine_centerpoints(output_folder_path, cache_folder_path, fragments_folder_path, consensus_swcs_folder_path, sample_folder_path) ;
