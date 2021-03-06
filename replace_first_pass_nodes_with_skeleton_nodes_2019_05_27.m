close all
clear

sample_date = '2019-05-27' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
%output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/%s/carver', sample_date) ;
intermediate_results_folder_path = fullfile(this_folder_path, sprintf('%s-first-pass-refined', sample_date)) ;
%fragments_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/prob0_swcs/frags', sample_date) ;
input_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;
sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
skeleton_folder_path =  sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;
%do_correct_fragments = true ;

% replace_machine_edges_with_skeleton_points(output_folder_path, ...
%                                            fragments_folder_path, ...
%                                            consensus_swcs_folder_path, ...
%                                            sample_folder_path, ...
%                                            skeleton_folder_path, ...
%                                            do_correct_fragments) ;

% Get the metadata for the rendered sample
[shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(sample_folder_path);

% % Load the fragments from disk                                                
% all_uncorrected_fragment_xyzs_in_erhan_coords = ...
%     compute_or_read_from_memo(output_folder_path, ...
%                               'all_uncorrected_fragment_xyzs_in_erhan_coords',  ...
%                               @()(read_all_fragment_centerpoints(fragments_folder_path))) ;

% % Check that the fragment nodes are on the grid as we expect
% %     all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
% %     all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
% %     fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
% %     assert( all( all( abs(fragment_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
% fraction_of_uncorrected_fragment_nodes_at_voxel_centers = ...
%     fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(all_uncorrected_fragment_xyzs_in_erhan_coords, spacing, origin)
% assert( fraction_of_uncorrected_fragment_nodes_at_voxel_centers == 1) ;
% 
% % Convert fragment coords to JaWS coords
% all_uncorrected_fragment_ijks_unrounded = (all_uncorrected_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
% all_uncorrected_fragment_ijks = round(all_uncorrected_fragment_ijks_unrounded) ;  % these should be one-based ijks

% Get the skeleton graph, either from .txt files or from mat                              
graph_file_path = fullfile(intermediate_results_folder_path, 'skeleton-as-graph.mat') ;
if exist(graph_file_path, 'file') ,
    load(graph_file_path, 'skeleton_graph', 'skeleton_ijks') ;
else        
    [skeleton_graph, skeleton_ijks] = load_skeleton_graph_from_txt_files(skeleton_folder_path, shape_xyz) ;  % these ijks are also one-based
    save(graph_file_path, 'skeleton_graph','skeleton_ijks', '-v7.3') ;
end

% % What fraction of uncorrected fragment points are also skeleton points
% is_uncorrected_fragment_a_skeleton_point = is_point_drawn_from_pool(all_uncorrected_fragment_ijks, skeleton_ijks, spacing) ;
% fraction_uncorrected_fragment_points_that_are_skeleton_points = mean(is_uncorrected_fragment_a_skeleton_point)
% 
% % Correct the fragment ijks, if requested
% if do_correct_fragments ,
%     all_fragment_ijks = all_uncorrected_fragment_ijks + 1 ;
% else
%     all_fragment_ijks = all_uncorrected_fragment_ijks ;        
% end   
% 
% % What fraction of corrected fragment points are also skeleton points
% is_fragment_a_skeleton_point = is_point_drawn_from_pool(all_fragment_ijks, skeleton_ijks, spacing) ;
% fraction_fragment_points_that_are_skeleton_points = mean(is_fragment_a_skeleton_point)
% assert( fraction_of_uncorrected_fragment_nodes_at_voxel_centers > 0.5) ;

% Convert one-based voxel indices to voxel center coords in JaWS coord system
% all_uncorrected_fragment_xyzs = jaws_origin + spacing .* (all_uncorrected_fragment_ijks-1) ;  % um, n x 3
% all_fragment_xyzs = jaws_origin + spacing .* (all_fragment_ijks-1) ;  % um, n x 3
skeleton_xyzs = jaws_origin + spacing .* (skeleton_ijks-1) ;  % um, n x 3
    % skeleton_ijks are 1-based, we want to map them to voxel centers
    % in JaWS coordinates

%     % Write the skeleton to json file, for Will & Jan
%     graph_as_json_file_path = fullfile(output_folder_path, 'skeleton-as-graph.json') ;
%     if ~exist(graph_as_json_file_path, 'file') ,
%         json_ready_skeleton = struct() ;
%         json_ready_skeleton.ijks = skeleton_ijks-1 ;  % convert to zero-based indexing
%         edges_with_third_col = table2array(skeleton_graph.Edges)-1 ;  % convert to zero-based indexing        
%         json_ready_skeleton.edges = edges_with_third_col(:,1:2) ;
%         json_ready_skeleton.jaws_origin = jaws_origin ;
%         json_ready_skeleton.spacing = spacing ;
%         json_ready_skeleton.origin = origin ;
%         json_ready_skeleton.shape_xyz = shape_xyz ;        
%         json_ready_skeleton.xyz_from_ijk_formula = 'xyzs = jaws_origin + spacing .* ijks' ;        
%         json_text = jsonencode(json_ready_skeleton) ;
%         dump_string_to_file(graph_as_json_file_path, json_text) ;
%     end    




% Load in the input neurons                          
input_neurons_and_names = compute_or_read_from_memo(intermediate_results_folder_path, ...
                                                    'input_neurons',  ...
                                                    @()(collect_first_pass_neurons(input_swcs_folder_path))) ;
input_neurons_as_swc_arrays =  input_neurons_and_names.input_neurons ;
input_neuron_names = input_neurons_and_names.input_neuron_names ;

% Check that those neurons are on-grid
% This currently fails for 2018-07-02 input neurons.  It's true
% that the vast majority of the nodes are on-grid, but some are not.
% Should investigate at some point
input_neuron_xyzs = xyzs_from_neuron_structs(input_neurons_as_swc_arrays) ;
%    assert( are_all_xyzs_at_voxel_centers(input_neuron_xyzs, spacing, jaws_origin) ) ;
fraction_of_input_nodes_at_voxel_centers = fraction_of_xyzs_at_voxel_centers_using_jaws_conventions(input_neuron_xyzs, spacing, jaws_origin)
%assert( fraction_of_input_nodes_at_voxel_centers == 1) ;

% % What fraction of the consensus nodes are fragment nodes for the uncorrected
% % fragments
% is_consensus_node_an_uncorrected_fragment_node = is_point_drawn_from_pool(consensus_neuron_xyzs, all_uncorrected_fragment_xyzs, spacing) ;
% fraction_of_consensus_nodes_that_are_uncorrected_fragment_nodes = mean(is_consensus_node_an_uncorrected_fragment_node)
% 
% % What fraction of the consensus nodes are fragment nodes
% is_consensus_node_a_fragment_node = is_point_drawn_from_pool(consensus_neuron_xyzs, all_fragment_xyzs, spacing) ;
% fraction_of_consensus_nodes_that_are_fragment_nodes = mean(is_consensus_node_a_fragment_node)

% %%
% % Label the consensus neuron nodes as being machine-generated or not
% consensus_neurons_with_machine_centerpoints_labelled = ...
%     compute_or_read_from_memo(output_folder_path, ...
%                               'consensus_neurons_with_machine_centerpoints_labelled',  ...
%                               @()(label_machine_centerpoints_in_neurons(consensus_neurons, consensus_neuron_names, all_fragment_xyzs, spacing))) ;                                  

%%                          
% Splice in skeleton nodes where possible
input_neurons_augmented_with_skeleton_nodes = ...
    compute_or_read_from_memo(intermediate_results_folder_path, ...
                              'input_neurons_augmented_with_skeleton_nodes',  ...
                              @()(replace_consensus_nodes_with_skeleton_nodes_in_neurons(input_neurons_as_swc_arrays, ...
                                                                                         input_neuron_names, ...
                                                                                         skeleton_graph, ...
                                                                                         skeleton_xyzs, ...
                                                                                         spacing))) ;                                  
                                                                                    
%%
% Write the swc's that don't exist already                          
input_neuron_count = length(input_neurons_as_swc_arrays) ;
swc_folder_path = fullfile(intermediate_results_folder_path, 'augmented-with-skeleton-nodes-as-swcs') ;
if consensus_neuron_count>0 && ~exist(swc_folder_path, 'file') ,
    mkdir(swc_folder_path) ;
end
for i = 1:input_neuron_count ,
    input_neuron_name = consensus_neuron_names{i} ;
    output_swc_file_name = sprintf('%s.swc', input_neuron_name) ;
    output_swc_file_path = fullfile(swc_folder_path, output_swc_file_name) ;
    if ~exist(output_swc_file_path, 'file') ,
        save_swc(output_swc_file_path, input_neurons_augmented_with_skeleton_nodes{i}, input_neuron_name) ;
    end
end
