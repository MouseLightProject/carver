function result = label_machine_centerpoints_in_neuron(neuron, kd_tree, machine_centerpoints, spacing)
    neuron_centerpoints = neuron(:,3:5) ;
    index_of_nearest_machine_point = knnsearch(kd_tree, neuron_centerpoints) ;
    nearest_machine_point = machine_centerpoints(index_of_nearest_machine_point, :) ;
    offset_to_nearest_machine_point = nearest_machine_point - neuron_centerpoints ;    
    offset_threshold = spacing/2 ;  % um 
    is_offset_below_threshold = bsxfun(@lt, offset_to_nearest_machine_point, offset_threshold) ;
    is_a_machine_point = all(is_offset_below_threshold, 2) ;  
    new_structure_identifier = 42 + is_a_machine_point ;
    result = neuron ;
    result(:,2) = new_structure_identifier ;
end