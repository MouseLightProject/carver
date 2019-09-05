function result = replace_machine_edges_with_skeleton_points_in_neurons(neurons, neuron_names, skeleton_graph, skeleton_xyzs, spacing)
    % neurons is a cell array, each element an centerpoint_count x 7 double array
    % Each row is an SWC-style record: id, type, x, y, z, something i forget,
    % and parent id.  The root has a parent id of -1, and parents come before
    % children.
    %
    % machine_centerpoints is a set of machine-annotated points, each row an x
    % y z triple
    
    skeleton_kd_tree  = KDTreeSearcher(skeleton_xyzs) ;        
    result = cell(size(neurons)) ;
    neuron_count = length(neurons) ;
    for i = 1:neuron_count ,
        neuron = neurons{i} ;
        neuron_name = neuron_names{i} ; 
        result{i} = ...
            replace_machine_edges_with_skeleton_points_in_neuron(neuron, neuron_name, skeleton_graph, skeleton_kd_tree, skeleton_xyzs, spacing) ;
    end
end
