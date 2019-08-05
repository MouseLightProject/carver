function result = label_machine_centerpoints_in_neurons(neurons, machine_centerpoints, spacing)
    % neurons is a cell array, each element an centerpoint_count x 7 double array
    % Each row is an SWC-style record: id, type, x, y, z, something i forget,
    % and parent id.  The root has a parent id of -1, and parents come before
    % children.
    %
    % machine_centerpoints is a set of machine-annotated points, each row an x
    % y z triple
    
    kd_tree  = KDTreeSearcher(machine_centerpoints) ;        
    result = cellfun(@(neuron)(label_machine_centerpoints_in_neuron(neuron, kd_tree, machine_centerpoints, spacing)), neurons, 'UniformOutput', false) ;
end
