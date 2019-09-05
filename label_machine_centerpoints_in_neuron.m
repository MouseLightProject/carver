function result = label_machine_centerpoints_in_neuron(neuron, neuron_name, kd_tree, machine_centerpoints, spacing)
%     if isequal(neuron_name, 'G-002') ,
%         fprintf('here!\n') ;
%     end
    neuron_centerpoints = neuron(:,3:5) ;
    index_of_nearest_machine_point = knnsearch(kd_tree, neuron_centerpoints) ;
    nearest_machine_point = machine_centerpoints(index_of_nearest_machine_point, :) ;
    offset_to_nearest_machine_point = nearest_machine_point - neuron_centerpoints ;    
    offset_threshold = spacing/16 ;  % um 
    is_offset_below_threshold = bsxfun(@lt, abs(offset_to_nearest_machine_point), offset_threshold) ;
    is_a_machine_point = all(is_offset_below_threshold, 2) ;  
    fprintf('Fraction machine points in neuron %s: %g\n', neuron_name, mean(is_a_machine_point)) ;
    new_structure_identifier = 42 + is_a_machine_point ;
    result = neuron ;
    result(:,2) = new_structure_identifier ;
end
