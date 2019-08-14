function result = are_all_neuron_nodes_on_jaws_grid(neuron, spacing, jaws_origin)
    % Returns true if all the nodes in the given neuron are on-grid, using
    % the conventions of Janelia Workstation circa 2019-08.
    % The neuron is represented in a zeven-column SWC-ready form.
    % spacing is 1 x 3, in um
    % jaws_origin is 1 x 3, in um
    neuron_xyzs = neuron(:, 3:5) ;  % um, xyz
    neuron_ijks_unrounded = (neuron_xyzs - jaws_origin) ./ spacing + 1 ;   % one-based indices
    neuron_ijks = round(neuron_ijks_unrounded) ;
    node_offsets = neuron_ijks_unrounded - neuron_ijks ;
    result =  all( all( node_offsets < spacing/16 ) ) ;  % the offsets should just be floating-point error
end
