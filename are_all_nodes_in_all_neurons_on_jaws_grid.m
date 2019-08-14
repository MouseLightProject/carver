function result = are_all_nodes_in_all_neurons_on_jaws_grid(neurons, spacing, jaws_origin)
    % Returns true if all the nodes in the given neurons are on-grid, using
    % the conventions of Janelia Workstation circa 2019-08.
    % Each neuron is represented in a seven-column SWC-ready form.
    % neurons is a cell array of these
    % spacing is 1 x 3, in um
    % jaws_origin is 1 x 3, in um
    
    result = all( cellfun(@(neuron)(are_all_neuron_nodes_on_jaws_grid(neuron, spacing, jaws_origin)), neurons) ) ;
end
