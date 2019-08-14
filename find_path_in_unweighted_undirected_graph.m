function [path, distance] = find_path_in_unweighted_undirected_graph(A, source_id, target_id)
    max_distance = 100 ;
    penultimate_front_set = [] ;
    last_front_set = source_id ;  % set at distance==0 from the source    
    for distance = 1:max_distance ,
        last_front_set_neighbors = neighbors(A, last_front_set) ;        
        front_set = setdiff(last_front_set_neighbors, union(last_front_set, penultimate_front_set)) ;  % these are at distance == distance
        if any(target_id==front_set) ,
            path = [] ;
            break
        end
        % If get here, set up for next iteration
        penultimate_front_set = last_front_set ;
        last_front_set = front_set ;
    end
end

function result = neighbors(A, set)
    rows_of_set = A(set,:) ;
    is_neighbor = any(rows_of_set, 1) ;
    result = find(is_neighbor) ;
end
