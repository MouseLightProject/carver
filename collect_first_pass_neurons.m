function result = collect_first_pass_neurons(input_swcs_folder_path)    
    [input_swc_file_paths, input_neuron_names] = collect_first_pass_swc_file_paths(input_swcs_folder_path) ;
    neurons = cellfun(@load_swc, input_swc_file_paths, 'UniformOutput', false) ;
    result = struct('neurons', {neurons}, ...
                    'neuron_names', {input_neuron_names}) ;        
end
