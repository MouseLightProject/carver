function [origin, spacing] = load_transform_txt(file_name)
    transform_as_raw_struct = configparser(file_name) ;
        % has fields, ox, oy, oz, sx, sy, sz, all in nm
        % also nl, the number of levels in the octree
    origin = [ transform_as_raw_struct.ox transform_as_raw_struct.oy transform_as_raw_struct.oz]/1e3 ;  % nm -> um
    spacing = [ transform_as_raw_struct.sx transform_as_raw_struct.sy transform_as_raw_struct.sz]/1e3/2^(transform_as_raw_struct.nl-1) ;  % nm -> um
        % determine the spacing of the voxels at the bottom level
    % Note that this is not *exactly* the origin used by Jaws as of this writing (2019-07-16).
    % But it's intended to be the location of the voxel *center* of the
    % 'lowest-index' voxel in the stack.
end
