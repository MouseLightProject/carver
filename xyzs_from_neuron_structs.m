function result = xyzs_from_neuron_structs(neuron_structs)
    xyzs_from_neuron_index = cellfun(@(neuron_struct)(xyzs_from_neuron_struct(neuron_struct)), neuron_structs, 'UniformOutput', false)  ;
    neuron_count = length(xyzs_from_neuron_index) ;
    result = zeros(0,3) ;
    
    for i = 1 : neuron_count ,
        result = vertcat(result, xyzs_from_neuron_index{i}) ; %#ok<AGROW>
    end
end
