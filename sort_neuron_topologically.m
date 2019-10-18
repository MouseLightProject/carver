function result = sort_neuron_topologically(swc_as_matrix)
    node_count = size(swc_as_matrix, 1) ;
    original_node_id_from_original_node_index = swc_as_matrix(:,1) ;
    if length(unique(original_node_id_from_original_node_index)) ~= node_count ,
        error('sort_neuron_topoloogically() requires the set of node ids to be [1, node_count]') ;
    end
    original_parent_node_id_from_original_node_index = swc_as_matrix(:,end) ;
    
    edges = [original_node_id_from_original_node_index original_parent_node_id_from_original_node_index] ;  % child, parent order
    edges(original_parent_node_id_from_original_node_index==-1, :) = [] ;
    
    dig_antirootward = digraph(edges(:,2), edges(:,1), ones(node_count-1,1), node_count) ;  % edges point away from the root
    new_order = toposort(dig_antirootward, 'Order', 'stable') ;
    %[new_order, dig_antirootward_sorted] = toposort(dig_antirootward, 'Order', 'stable') ;
    swc_as_matrix_in_new_order = swc_as_matrix(new_order,:) ;
    original_node_id_from_new_node_index = swc_as_matrix_in_new_order(:,1) ;
    original_parent_node_id_from_new_node_index = swc_as_matrix_in_new_order(:,end) ;  % has a -1 for the root node
    is_root_node_from_new_node_index = (original_parent_node_id_from_new_node_index==-1) ;
    original_parent_node_id_from_new_node_index(is_root_node_from_new_node_index) = 1 ;  % patch so we can use it for indexing
    new_node_index_from_original_node_id = invert_map_array(original_node_id_from_new_node_index) ;
    new_node_id_from_new_node_index = new_node_index_from_original_node_id(original_node_id_from_new_node_index) ;  % should be identity permuation   
    assert( isequal(new_node_id_from_new_node_index, (1:node_count)') ) ;
    new_parent_node_id_from_new_node_index = new_node_index_from_original_node_id(original_parent_node_id_from_new_node_index) ;    
    new_parent_node_id_from_new_node_index(is_root_node_from_new_node_index) = -1 ;  % put -1 for the parent of the root node
    
    %new_parent_node_id_from_node_index_check_as_cell_array = ...
    %    arrayfun(@(node_index)(dig_antirootward_sorted.predecessors(node_index)), (1:node_count)', ...
    %             'UniformOutput', false) ;
    %new_parent_node_id_from_node_index_check = cellfun(@scalar_from_array_with_empty_encoded, new_parent_node_id_from_node_index_check_as_cell_array) ;
    %assert( isequal(new_parent_node_id_from_new_node_index, new_parent_node_id_from_node_index_check) ) ;    
    
    result = swc_as_matrix_in_new_order ;
    result(:,1) = new_node_id_from_new_node_index ;
    result(:,end) = new_parent_node_id_from_new_node_index ;
end



% function result = scalar_from_array_with_empty_encoded(array)
%     % This is used to convert the output of digraph::predecessors() to
%     % a parent array like is used in swc, using the swc convention that the
%     % parent of the root is -1.
%     %
%     if isempty(array)
%         result = -1 ;
%     elseif isscalar(array)
%         result = array ;
%     else
%         error('Inputs to scalar_from_array_with_empty_encoded() must be empty or scalar') ;
%     end    
% end
