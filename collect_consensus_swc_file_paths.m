function [consensus_swc_paths, neuron_names_with_consensus_swc_paths]  = collect_consensus_swc_file_paths(consensus_swcs_folder_path)    
    % E.g. one consensus .swc file path within the consensus_swcs_folder_path
    %
    % /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-10-01
    %
    % is
    % 
    % /groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-10-01/G-014/consensus/2018-10-01_G-014_consensus.swc
    
    date_as_string = leaf_file_name_from_path(consensus_swcs_folder_path) ;
    neuron_folder_path_template = fullfile(consensus_swcs_folder_path, 'G-*') ;
    neuron_names = simple_dir(neuron_folder_path_template) ;
    
    function result = consensus_swc_path_from_neuron_name(neuron_name)
        lower_case_file_name = sprintf('%s_%s_consensus.swc', date_as_string, neuron_name) ;
        lower_case_file_path = fullfile(consensus_swcs_folder_path, neuron_name, 'consensus', lower_case_file_name) ;
        if exist(lower_case_file_path, 'file') ,
            result = lower_case_file_path ;
        else
            upper_case_file_name = sprintf('%s_%s_Consensus.swc', date_as_string, neuron_name) ;
            upper_case_file_path = fullfile(consensus_swcs_folder_path, neuron_name, 'Consensus', upper_case_file_name) ;
            if exist(upper_case_file_path, 'file') ,
                result = upper_case_file_path ;
            else
                result = '' ;
            end
        end
    end

    consensus_swc_paths_but_with_empties = cellfun(@consensus_swc_path_from_neuron_name, neuron_names, 'UniformOutput', false) ;
    has_consensus_swc = cellfun(@is_nonempty, consensus_swc_paths_but_with_empties) ;
    consensus_swc_paths = consensus_swc_paths_but_with_empties(has_consensus_swc) ;
    neuron_names_with_consensus_swc_paths = neuron_names(has_consensus_swc) ;
end
