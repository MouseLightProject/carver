function result = find_path_greedily(A, xyzs, source_id, target_id)
    max_path_nodes = 1000 ;
    target_xyz = xyzs(target_id, :) ;
    did_find_path = false ;
    path_so_far = zeros(1, max_path_nodes) ;
    path_so_far(1) = source_id ;
    for path_node_index = 2:max_path_nodes ,
        current_node_id = path_so_far(path_node_index-1) ;
        neighbor_ids = get_neighbors(A, current_node_id) ;        
        neighbor_xyzs = xyzs(neighbor_ids,:) ;
        distance_to_target_per_neighbor = vecnorm(neighbor_xyzs - target_xyz, 2, 2) ;
        [~, j_best] = min(distance_to_target_per_neighbor) ;
        this_path_node_id = neighbor_ids(j_best) ;
        path_so_far(path_node_index) = this_path_node_id ;
        if this_path_node_id == target_id ,
            did_find_path = true ;
            path_node_count = path_node_index ;
            break
        end
    end
    if did_find_path ,
        result = path_so_far(1:path_node_count) ;
    else
        error('Unable to find path to target with %d nodes or fewer', max_path_nodes) ;
    end
end

function result = get_neighbors(A, set)
    rows_of_set = A(set, :) ;
    is_neighbor = any(rows_of_set, 1) ;
    result = find(is_neighbor) ;
end
