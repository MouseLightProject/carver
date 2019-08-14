function [did_succeed, path] = find_path_via_a_star_in_graph_proper(G, xyzs, source_id, target_id)
    source_xyz = xyzs(source_id, :) ;
    target_xyz = xyzs(target_id, :) ;
    distance_from_source_to_target = vecnorm(source_xyz-target_xyz, 2, 2) ;
    node_id_from_visit_index = source_id ;
    is_in_open_set_from_visit_index = true ;
    is_in_closed_set_from_visit_index = false ;
    predecessor_visit_index_from_visit_index = nan ;
    cheapest_path_known_cost_from_visit_index = 0 ;  % gScore, cheapest path known from source to this node
    heuristic_cost_to_target_from_visit_index = distance_from_source_to_target ;   % fScore, estimated total path to target 
    
    while any(is_in_open_set_from_visit_index) ,
        % The current node is the one in the open set with cheapest known
        % path from the source
        modified_heuristic_cost_to_target_from_visit_index = heuristic_cost_to_target_from_visit_index ;
        modified_heuristic_cost_to_target_from_visit_index(~is_in_open_set_from_visit_index) = inf ;    
        [heuristic_cost_to_target_from_current_node, current_visit_index] = min(modified_heuristic_cost_to_target_from_visit_index) ;
        if heuristic_cost_to_target_from_current_node > 10*distance_from_source_to_target ,
            % if the heuristic cost for the best-bet node right now is much
            % bigger than the straight-line distance from source to target,
            % it likely means there is no path in the skeleton that
            % approximates the source-to-target path.  This sometimes
            % happens, maybe b/c tracers just manually connect two
            % machine-generated points.  In this case, just punt and return
            % failure.
            fprintf('Failed to find a skeleton path!\n') ;
            did_succeed = false ;            
            path = zeros(1,0) ;
            return
        end
        current_node_id = node_id_from_visit_index(current_visit_index) ;
        
        if current_node_id == target_id , 
            did_succeed = true ;     
            path = reconstruct_path(predecessor_visit_index_from_visit_index, node_id_from_visit_index, current_visit_index) ;
            return
        end
                    
        is_in_open_set_from_visit_index(current_visit_index) = false ;
        is_in_closed_set_from_visit_index(current_visit_index) = true ;
        
        node_id_from_neighbor_index = (G.neighbors(current_node_id))' ;
        is_same_from_neighbor_index_and_visit_index = (node_id_from_visit_index == node_id_from_neighbor_index') ;
        
        % Create and populate visit_index_from_neighbor_index
        visit_index_from_neighbor_index = zeros(size(node_id_from_neighbor_index)) ;
        for neighbor_index = 1:length(node_id_from_neighbor_index) ,
            for visit_index = 1:length(node_id_from_visit_index) ,
                if is_same_from_neighbor_index_and_visit_index(neighbor_index, visit_index) ,
                    visit_index_from_neighbor_index(neighbor_index) = visit_index ;
                    break  % out of inner loop only
                end
            end
        end
                
        for neighbor_index = 1:length(node_id_from_neighbor_index) ,
            neighbor_node_id = node_id_from_neighbor_index(neighbor_index) ;
            neighbor_visit_index = visit_index_from_neighbor_index(neighbor_index) ;
            
            % If this node has not been visited before, add an entry for it
            % to the arrays that keep track of visited nodes
            if neighbor_visit_index == 0 ,
                % First time visiting this node
                neighbor_visit_index = length(node_id_from_visit_index) + 1 ;                
                node_id_from_visit_index(neighbor_visit_index) = neighbor_node_id ;
                is_in_open_set_from_visit_index(neighbor_visit_index) = true ;
                is_in_closed_set_from_visit_index(neighbor_visit_index) = false ;
                predecessor_visit_index_from_visit_index(neighbor_visit_index) = nan ;
                cheapest_path_known_cost_from_visit_index(neighbor_visit_index) = inf ;
                heuristic_cost_to_target_from_visit_index(neighbor_visit_index) = inf ;   % fScore, estimated total path to target 
            end
            
            is_in_closed_set = is_in_closed_set_from_visit_index(neighbor_visit_index) ;
            if is_in_closed_set ,
                continue 
            end
            distance_from_current_to_neighbor =  vecnorm(xyzs(neighbor_node_id, :)-xyzs(current_node_id, :), 2, 2) ;  
            tentative_cheapest_path_known_cost = cheapest_path_known_cost_from_visit_index(current_visit_index) + distance_from_current_to_neighbor ;
            if tentative_cheapest_path_known_cost < cheapest_path_known_cost_from_visit_index(neighbor_visit_index) ,
                % This path is best so far, so record it
                predecessor_visit_index_from_visit_index(neighbor_visit_index) = current_visit_index ;
                cheapest_path_known_cost_from_visit_index(neighbor_visit_index) = tentative_cheapest_path_known_cost ;
                distance_from_neighbor_to_target = vecnorm(xyzs(neighbor_node_id, :)-target_xyz, 2, 2) ;  
                heuristic_cost_to_target_from_visit_index(neighbor_visit_index) = tentative_cheapest_path_known_cost + distance_from_neighbor_to_target ;
            end
        end
    end
    
    % If get here, open set is empty, but target not reached
    fprintf('Failed to find a skeleton path!\n') ;
    did_succeed = false ;            
    path = zeros(1,0) ;
end



function result = reconstruct_path(predecessor_visit_index_from_visit_index, node_id_from_visit_index, target_visit_index)
    visit_index_from_path_index = target_visit_index ;    
    predecessor_visit_index = predecessor_visit_index_from_visit_index(target_visit_index) ;
    while isfinite(predecessor_visit_index) ,
        visit_index_from_path_index = horzcat(predecessor_visit_index, visit_index_from_path_index) ; %#ok<AGROW>
        % Get next predecessor
        predecessor_visit_index = predecessor_visit_index_from_visit_index(predecessor_visit_index) ;
    end
    result = node_id_from_visit_index(visit_index_from_path_index) ;
end


%     while openSet is not empty
%         current := the node in openSet having the lowest fScore[] value
%         if current = goal
%             return reconstruct_path(cameFrom, current)
% 
%         openSet.Remove(current)
%         closedSet.Add(current)
%         for each neighbor of current
%             if neighbor in closedSet 
%                 continue
%             % d(current,neighbor) is the weight of the edge from current to neighbor
%             % tentative_gScore is the distance from start to the neighbor through current
%             tentative_gScore := gScore[current] + d(current, neighbor)
%             if neighbor not in openSet
%                 openSet.add(neighbor)
%             if tentative_gScore < gScore[neighbor]
%                 % This path to neighbor is better than any previous one. Record it!
%                 cameFrom[neighbor] := current
%                 gScore[neighbor] := tentative_gScore
%                 fScore[neighbor] := gScore[neighbor] + h(neighbor)
% 
%     % Open set is empty but goal was never reached
%     return failure
    