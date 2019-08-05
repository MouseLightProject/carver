function result = collect_consensus_neurons(consensus_swcs_folder_path)    
    [consensus_swc_file_paths, consensus_neuron_names] = collect_consensus_swc_file_paths(consensus_swcs_folder_path) ;
    consensus_neurons = cellfun(@load_swc, consensus_swc_file_paths, 'UniformOutput', false) ;
    result = struct('consensus_neurons', {consensus_neurons}, ...
                    'consensus_neuron_names', {consensus_neuron_names}) ;        
end
