function result = replace_consensus_nodes_with_skeleton_nodes_in_neurons(consensus_neurons_as_swc_arrays, ...
                                                                         consensus_neuron_names, ...
                                                                         skeleton_graph, ...
                                                                         skeleton_xyzs, ...
                                                                         spacing)
    % neurons is a cell array, each element an centerpoint_count x 7 double array
    % Each row is an SWC-style record: id, type, x, y, z, something i forget,
    % and parent id.  The root has a parent id of -1, and parents come before
    % children.
    %
    % machine_centerpoints is a set of machine-annotated points, each row an x
    % y z triple
    
    skeleton_kd_tree  = KDTreeSearcher(skeleton_xyzs) ;        
    result = cell(size(consensus_neurons_as_swc_arrays)) ;
    neuron_count = length(consensus_neurons_as_swc_arrays) ;
    for i = 1:neuron_count ,
        neuron = consensus_neurons_as_swc_arrays{i} ;
        neuron_name = consensus_neuron_names{i} ; 
        result{i} = ...
            replace_consensus_nodes_with_skeleton_nodes_in_neuron(neuron, neuron_name, skeleton_graph, skeleton_kd_tree, skeleton_xyzs, spacing) ;
    end
end
