function result = fraction_of_xyzs_at_voxel_centers(xyzs, spacing, origin)
    % Returns the fraction of points at voxel centers.  Assumes that origin is
    % the location of *center* of the voxel with lowest-possible indices.
    % This is the convention used by the  Janelia Workstation circa 2019-08.
    %
    % xyzs should be node_count x 3, in um, and in xyz order
    % spacing is 1 x 3, in um
    % origin is 1 x 3, in um
    ijks_unrounded = (xyzs - origin) ./ spacing + 1 ;   % one-based indices
    ijks = round(ijks_unrounded) ;
    node_offsets = ijks_unrounded - ijks ;
    absolute_node_offsets = abs(node_offsets) ;
    is_at_voxel_center_from_neuron_index = all( absolute_node_offsets < 0.001, 2 ) ;  % the offsets should just be floating-point error
    result = mean( is_at_voxel_center_from_neuron_index ) ;
end
