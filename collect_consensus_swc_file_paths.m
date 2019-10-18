function [consensus_swc_paths, neuron_names_with_consensus_swc_paths]  = collect_consensus_swc_file_paths(consensus_swcs_root_folder_path)    
    % E.g. one consensus .swc file path within the consensus_swcs_folder_path
    %
    % /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-10-01
    %
    % is
    % 
    % /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-10-01/G-014/consensus/2018-10-01_G-014_consensus.swc
    
    %date_as_string = leaf_file_name_from_path(consensus_swcs_root_folder_path) ;
    neuron_folder_path_template = fullfile(consensus_swcs_root_folder_path, 'G-*') ;
    neuron_names = simple_dir(neuron_folder_path_template) ;
    
    function result = consensus_swc_path_from_neuron_name(neuron_name)
        result = '' ;  % fallback if early return
        lower_case_folder_path = fullfile(consensus_swcs_root_folder_path, neuron_name, 'consensus') ;
        upper_case_folder_path = fullfile(consensus_swcs_root_folder_path, neuron_name, 'Consensus') ;
        possible_consensus_folder_paths = {lower_case_folder_path upper_case_folder_path} ;
        does_folder_path_exist = cellfun(@(path)(logical(exist(path, 'file'))), possible_consensus_folder_paths) ;
        folder_path_count = sum(does_folder_path_exist) ;
        if folder_path_count == 0 ,
            fprintf('Warning: Can''t find a consensus folder for neuron %s\n', neuron_name) ;
            return
        elseif folder_path_count > 1 ,
            fprintf('Warning: Both "consensus" and "Consensus" folders exist for neuron %s.  Don''t know which one to use, so skipping this neuron.\n', neuron_name) ;
            return
        end
        consensus_folder_path = possible_consensus_folder_paths{does_folder_path_exist} ;
        
        swc_file_path_path_template = fullfile(consensus_folder_path, '*onsensus.swc') ;
        consensus_swc_file_names = simple_dir(swc_file_path_path_template) ;
        consensus_swc_file_name_count = length(consensus_swc_file_names) ;
        if consensus_swc_file_name_count == 0 ,
            fprintf('Warning: Can''t find a consensus .swc file in folder %s\n', consensus_folder_path) ;
            return
        elseif folder_path_count > 1 ,
            fprintf('Warning: More than one consensus .swc file in folder %s.  Don''t know which one to use, so skipping this neuron.\n', consensus_folder_path) ;
            return
        end        
        consensus_swc_file_name = consensus_swc_file_names{1} ;
        result = fullfile(consensus_folder_path, consensus_swc_file_name) ;
    end

    consensus_swc_paths_but_with_empties = cellfun(@consensus_swc_path_from_neuron_name, neuron_names, 'UniformOutput', false) ;
    has_consensus_swc = cellfun(@is_nonempty, consensus_swc_paths_but_with_empties) ;
    if ~all(has_consensus_swc) ,
        fprintf('Some consensus neurons had a folder, but seemingly no consensus .swc inside:\n') ;
        missing_neuron_names = neuron_names(~has_consensus_swc)   %#ok<NOPRT,NASGU>
    end
    consensus_swc_paths = consensus_swc_paths_but_with_empties(has_consensus_swc) ;
    neuron_names_with_consensus_swc_paths = neuron_names(has_consensus_swc) ;
end
