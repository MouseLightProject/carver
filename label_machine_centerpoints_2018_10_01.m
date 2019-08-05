output_folder_path = '/home/taylora/cache/gt/2018-10-01' ;
cache_folder_path = '/home/taylora/cache/gt/2018-10-01' ;
fragments_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2018-10-01/prob0_swcs/frags' ;
consensus_swcs_folder_path = '/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-10-01' ;
sample_folder_path = '/nrs/mouselight/SAMPLES/2018-10-01' ;

consensus_neurons_with_machine_centerpoints_labelled = ...
    label_machine_centerpoints(output_folder_path, cache_folder_path, fragments_folder_path, consensus_swcs_folder_path, sample_folder_path) ;
