function result = replace_machine_edges_with_skeleton_points_in_neuron(consensus_neuron_as_swc_array, ...
                                                                       consensus_neuron_name, ...
                                                                       skeleton_graph, ...
                                                                       skeleton_kd_tree, ...
                                                                       skeleton_xyzs, ...
                                                                       spacing)
    xyz_from_node_id = consensus_neuron_as_swc_array(:,3:5) ;
    is_a_fragment_point_from_node_id = (consensus_neuron_as_swc_array(:,2)==43) ;
    parent_node_id_from_node_id = consensus_neuron_as_swc_array(:, end) ;
    if ~isempty(parent_node_id_from_node_id) ,
        parent_node_id_from_node_id(1) = 1 ;  % first element is the root, but convenient to have 'parent' be a valid node id
    end
        
    index_of_nearest_skeleton_point = knnsearch(skeleton_kd_tree, xyz_from_node_id) ;
    nearest_skeleton_xyz = skeleton_xyzs(index_of_nearest_skeleton_point, :) ;
    offset_to_nearest_skeleton_point = nearest_skeleton_xyz - xyz_from_node_id ;    
    absolute_offset_to_nearest_skeleton_point = abs(offset_to_nearest_skeleton_point) ;
    % offset_threshold = spacing/16 ;  % um 
    offset_threshold = 1.5 * spacing;  % um, 1x3 
    is_offset_below_threshold = (absolute_offset_to_nearest_skeleton_point < offset_threshold) ;
    is_a_periskeleton_point_from_node_id = all(is_offset_below_threshold, 2) ;
    
    is_parent_a_fragment_point_from_node_id = is_a_fragment_point_from_node_id(parent_node_id_from_node_id) ;
    is_parent_a_fragment_point_from_node_id(1) = false ;  % patch root node
    is_self_and_parent_a_fragment_point_from_node_id = ... 
        is_a_fragment_point_from_node_id & is_parent_a_fragment_point_from_node_id ;

    is_parent_a_periskeleton_point_from_node_id = is_a_periskeleton_point_from_node_id(parent_node_id_from_node_id) ;
    is_parent_a_periskeleton_point_from_node_id(1) = false ;  % patch root node
    is_self_and_parent_a_periskeleton_point_from_node_id = ... 
        is_a_periskeleton_point_from_node_id & is_parent_a_periskeleton_point_from_node_id ;
    
    is_self_and_parent_both_fragment_and_periskeleton_from_node_id = ...
        is_self_and_parent_a_fragment_point_from_node_id & is_self_and_parent_a_periskeleton_point_from_node_id ;
    
%     if isequal(neuron_name, 'G-002') ,
%         fprintf('here!\n') ;
%     end
    node_ids_to_replace = find(is_self_and_parent_both_fragment_and_periskeleton_from_node_id) ;
    %A_skeleton_graph = skeleton_graph.adjacency ;
    node_ids_to_replace_count = length(node_ids_to_replace) ;
    is_a_skeleton_point_from_node_id = false(size(is_a_fragment_point_from_node_id)) ;
    progress_bar('reset') ;
    for i = 1 : node_ids_to_replace_count ,
        node_id = node_ids_to_replace(i) ;
        parent_node_id = parent_node_id_from_node_id(node_id) ;
        skeleton_node_id = index_of_nearest_skeleton_point(node_id) ;
        parent_skeleton_node_id = index_of_nearest_skeleton_point(parent_node_id) ;    
%         tic_id = tic() ;
%         %[~, skeleton_path] = graphshortestpath(A_skeleton_graph, parent_skeleton_node_id, skeleton_node_id, 'directed', false, 'Method', 'BFS') ;
%         %skeleton_path = find_path_greedily_in_graph_proper(skeleton_graph, skeleton_xyzs , parent_skeleton_node_id, skeleton_node_id) ;
%         skeleton_path = skeleton_graph.shortestpath(parent_skeleton_node_id, skeleton_node_id, 'Method', 'unweighted') ;
%         elapsed_time = toc(tic_id) ;
%         fprintf('Elapsed time for BFS search: %g seconds\n', elapsed_time) ;
        %tic_id = tic() ;
        [did_find_skeleton_path, skeleton_path] = ...
            find_path_via_a_star_in_graph_proper(skeleton_graph, skeleton_xyzs, parent_skeleton_node_id, skeleton_node_id) ;
        %elapsed_time = toc(tic_id) ;
        %fprintf('Elapsed time for A* search: %g seconds\n', elapsed_time) ;
        %assert( isequal(skeleton_path, skeleton_path_2) ) ;
        
        if did_find_skeleton_path ,
            xyz_from_node_id(node_id,:) = skeleton_xyzs(skeleton_node_id, :) ;  
                % substitute the skeleton node for the node (although the skeleton node will
                % often be identical)
            is_a_skeleton_point_from_node_id(node_id) = true ;
            xyz_from_node_id(parent_node_id,:) = skeleton_xyzs(parent_skeleton_node_id, :) ;
                % substitute the parent skeleton node for the parent node (although the skeleton node will
                % often be identical)            
            is_a_skeleton_point_from_node_id(parent_node_id) = true ;
            trimmed_skeleton_path = skeleton_path(2:end-1) ;  % trim off the start and end nodes, since we already handled those
            new_node_count = length(trimmed_skeleton_path) ;
            if new_node_count>0 ,
                new_node_xyzs = skeleton_xyzs(trimmed_skeleton_path, :) ;
                % Add the new nodes at end the end of the node list
                node_count = length(parent_node_id_from_node_id) ;
                new_node_ids = node_count + (1:new_node_count)' ;
                new_parent_node_id_from_node_id = [parent_node_id ; new_node_ids(1:end-1)] ;
                xyz_from_node_id = [ xyz_from_node_id ; new_node_xyzs] ;  %#ok<AGROW>
                parent_node_id_from_node_id = [parent_node_id_from_node_id ; new_parent_node_id_from_node_id ] ;  %#ok<AGROW>
                parent_node_id_from_node_id(node_id) = new_node_ids(end) ;
                is_a_skeleton_point_from_node_id = [ is_a_skeleton_point_from_node_id ; true(new_node_count, 1) ] ;  %#ok<AGROW>
            end
        end

        progress_bar(i, node_ids_to_replace_count) ;
    end
    
    % Put stuff together into a neuron_swc_array
    result_node_count = length(parent_node_id_from_node_id) ;
    parent_node_id_from_node_id(1) = -1 ;  % use .swc convention for parent of the root node
    unsorted_result = [ (1:result_node_count)' 52+is_a_skeleton_point_from_node_id xyz_from_node_id ones(result_node_count,1) parent_node_id_from_node_id ] ;
    
    % Do a topological sort on that, so parents always come before children
    result = sort_neuron_topologically(unsorted_result) ;
end
