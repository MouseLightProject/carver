sample_date = '2018-08-01' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
output_folder_path = fullfile(this_folder_path, sprintf('%s-consensus-swcs', sample_date)) ;
consensus_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;

% Copy the consensus neurons                          
snapshot_consensus_neurons(consensus_swcs_folder_path, output_folder_path) ;
