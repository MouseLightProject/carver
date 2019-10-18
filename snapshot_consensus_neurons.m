function snapshot_consensus_neurons(consensus_swcs_folder_path, output_folder_path)    
    [consensus_swc_file_paths, consensus_neuron_names] = collect_consensus_swc_file_paths(consensus_swcs_folder_path) ;
    swc_count = length(consensus_swc_file_paths) ;
    if ~exist(output_folder_path, 'file') ,
        mkdir(output_folder_path) ;
    end
    for i = 1 : swc_count ,
        consensus_swc_file_path = consensus_swc_file_paths{i} ;
        consensus_swc_file_name = file_name_from_path(consensus_swc_file_path) ;
        output_swc_file_path = fullfile(output_folder_path, consensus_swc_file_name) ;
        if ~exist(output_swc_file_path, 'file') ,
            copyfile(consensus_swc_file_path, output_swc_file_path) ;
        end
    end        
end
