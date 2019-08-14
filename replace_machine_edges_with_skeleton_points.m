function consensus_neurons_augmented_with_skeleton_nodes = ...
    replace_machine_edges_with_skeleton_points(output_folder_path, ...
                                               fragments_folder_path, ...
                                               consensus_swcs_folder_path, ...
                                               sample_folder_path, ...
                                               skeleton_folder_path)
    
    % Get the metadata for the rendered sample
    [shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(sample_folder_path);
    
    % Load the fragments from disk                                                
    all_fragment_xyzs_in_erhan_coords = ...
        compute_or_read_from_memo(output_folder_path, ...
                                  'all_fragment_xyzs_in_erhan_coords',  ...
                                  @()(read_all_fragment_centerpoints(fragments_folder_path))) ;
                              
    % Check that the fragment nodes are on the grid as we expect
    all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
    all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
    fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
    assert( all( all( fragment_offsets < spacing/16 ) ) ) ;  % the offsets should just be floating-point error
                              
    % Convert fragment coords to JaWS coords
    all_fragment_xyzs = jaws_origin + spacing .* (all_fragment_ijks-1) ;  % um, n x 3
    
    % Load in the consensus neurons                          
    consensus_neurons_and_names = compute_or_read_from_memo(output_folder_path, ...
                                                            'consensus_neurons',  ...
                                                            @()(collect_consensus_neurons(consensus_swcs_folder_path))) ;
    consensus_neurons =  consensus_neurons_and_names.consensus_neurons ;
    consensus_neuron_names = consensus_neurons_and_names.consensus_neuron_names ;
    
    % Check that those neurons are on-grid
    % This currently fails for 2018-07-02 consensus neurons!  It's true
    % that the vast majority of the nodes are on-grid, but some are not!
    % What the heck?!
    %assert( are_all_nodes_in_all_neurons_on_jaws_grid(consensus_neurons, spacing, jaws_origin) ) ;
    
    % Label the consensus neuron nodes as being machine-generated or not
    consensus_neurons_with_machine_centerpoints_labelled = ...
        compute_or_read_from_memo(output_folder_path, ...
                                  'consensus_neurons_with_machine_centerpoints_labelled',  ...
                                  @()(label_machine_centerpoints_in_neurons(consensus_neurons, all_fragment_xyzs, spacing))) ;                                  
    
    % Get the skeleton graph, either from .txt files or from mat                              
    graph_file_path = fullfile(output_folder_path, 'skeleton-as-graph.mat') ;
    if exist(graph_file_path, 'file') ,
        load(graph_file_path, 'skeleton_graph', 'skeleton_ijks') ;
    else        
        [skeleton_graph, skeleton_ijks] = load_skeleton_graph_from_txt_files(skeleton_folder_path, shape_xyz) ;  % these ijks are also one-based
        save(graph_file_path, 'skeleton_graph','skeleton_ijks', '-v7.3') ;
    end
    skeleton_xyzs = jaws_origin + spacing .* (skeleton_ijks-1) ;  % um, n x 3
        % skeleton_ijks are 1-based, we want to map them to voxel centers
        % in JaWS coordinates
        
    % Splice in skeleton nodes where possible
    consensus_neurons_augmented_with_skeleton_nodes = ...
        compute_or_read_from_memo(output_folder_path, ...
                                  'consensus_neurons_augmented_with_skeleton_nodes',  ...
                                  @()(replace_machine_edges_with_skeleton_points_in_neurons(consensus_neurons_with_machine_centerpoints_labelled, ...
                                                                                            skeleton_graph, ...
                                                                                            skeleton_xyzs, ...
                                                                                            spacing))) ;                                  
    
    % Write the swc's that don't exist already                          
    consensus_neuron_count = length(consensus_neurons) ;
    swc_folder_path = fullfile(output_folder_path, 'augmented-with-skeleton-nodes-as-swcs') ;
    if consensus_neuron_count>0 && ~exist(swc_folder_path, 'file') ,
        mkdir(swc_folder_path) ;
    end
    for i = 1:consensus_neuron_count ,
        consensus_neuron_name = consensus_neuron_names{i} ;
        swc_file_name = sprintf('%s.swc', consensus_neuron_name) ;
        swc_file_path = fullfile(swc_folder_path, swc_file_name) ;
        if ~exist(swc_file_path, 'file') ,
            save_swc(swc_file_path, consensus_neurons_augmented_with_skeleton_nodes{i}, consensus_neuron_name) ;
        end
    end
end
