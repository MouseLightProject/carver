function [swc_paths, extended_neuron_names]  = collect_first_pass_swc_file_paths(input_swcs_folder_path)    
    % Collects .swc file names within a e.g. G-013 folders within
    % input_swcs_folder_path.
    neuron_folder_path_template = fullfile(input_swcs_folder_path, 'G-*') ;
    neuron_names = simple_dir(neuron_folder_path_template) ;
    
    function result = first_pass_swc_paths_from_neuron_name(neuron_name)
        neuron_folder_path = fullfile(input_swcs_folder_path, neuron_name) ;
        swc_file_path_template = fullfile(neuron_folder_path, '*.swc') ;
        swc_file_names = simple_dir(swc_file_path_template) ;
        result = cellfun(@(swc_file_name)(fullfile(neuron_folder_path, swc_file_name)), ...
                         swc_file_names, ...
                         'UniformOutput', false) ;
    end

    list_of_list_of_input_swc_paths = cellfun(@first_pass_swc_paths_from_neuron_name, neuron_names, 'UniformOutput', false) ;

    % flatten the list of lists
    swc_paths = cell(1,0) ;
    n = length(list_of_list_of_input_swc_paths) ;
    for i = 1:n ,
        swc_paths_for_this_neuron = list_of_list_of_input_swc_paths{i} ;
        swc_paths = horzcat(swc_paths, swc_paths_for_this_neuron) ;  %#ok<AGROW>
    end    
    
    extended_neuron_names = cellfun(@leaf_base_name_from_path, swc_paths, 'UniformOutput', false) ;    
end



function result = leaf_base_name_from_path(path)
    [~, result] = fileparts(path) ;
end
