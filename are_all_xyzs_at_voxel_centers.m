function result = are_all_xyzs_at_voxel_centers(xyzs, spacing, jaws_origin)
    % Returns true if all the points given are on-grid, using
    % the conventions of Janelia Workstation circa 2019-08.
    % xyzs should be node_count x 3, in um, and in xyz order
    % spacing is 1 x 3, in um
    % jaws_origin is 1 x 3, in um
    ijks_unrounded = (xyzs - jaws_origin) ./ spacing + 1 ;   % one-based indices
    ijks = round(ijks_unrounded) ;
    node_offsets = ijks_unrounded - ijks ;
    absolute_node_offsets = abs(node_offsets) ;
    result =  all( all( absolute_node_offsets < 1/16 ) ) ;  % the offsets should just be floating-point error
end
