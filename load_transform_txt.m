function [origin, spacing, octree_level_count] = load_transform_txt(file_name)
    transform_as_raw_struct = configparser(file_name) ;
        % has fields, ox, oy, oz, sx, sy, sz, all in nm
        % also nl, the number of levels in the octree
    octree_level_count = transform_as_raw_struct.nl ;
    octree_level_steps = octree_level_count - 1 ;
    origin = [ transform_as_raw_struct.ox transform_as_raw_struct.oy transform_as_raw_struct.oz]/1e3 ;  % nm -> um
    spacing = [ transform_as_raw_struct.sx transform_as_raw_struct.sy transform_as_raw_struct.sz]/1e3/2^octree_level_steps ;  % nm -> um
        % determine the spacing of the voxels at the bottom level
    % Note that this is not *exactly* the origin used by Jaws as of this writing (2019-07-16).
    % But it's intended to be the location of the voxel *center* of the
    % 'lowest-index' voxel in the stack.
    
end
