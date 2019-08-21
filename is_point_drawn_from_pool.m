function result = is_point_drawn_from_pool(test_points, pool_points, spacing)
    % points is test_point_count x 3
    % pool is pool_point_count x 3
    % spacing is 1 x 3, represents the voxel spacing that points and pool
    % should be on
    kd_tree  = KDTreeSearcher(pool_points) ;        
    index_of_nearest_pool_point = knnsearch(kd_tree, test_points) ;
    nearest_pool_point = pool_points(index_of_nearest_pool_point, :) ;
    offset_to_nearest_pool_point = nearest_pool_point - test_points ;    
    normed_offset_to_nearest_pool_point = offset_to_nearest_pool_point ./ spacing ;
    absolute_normed_offset_to_nearest_pool_point = abs(normed_offset_to_nearest_pool_point) ;
    max_absolute_normed_offset_to_nearest_pool_point = max(absolute_normed_offset_to_nearest_pool_point, [], 2) ;
    normed_offset_threshold = 1/2 ;
    result = (max_absolute_normed_offset_to_nearest_pool_point < normed_offset_threshold) ;
end
