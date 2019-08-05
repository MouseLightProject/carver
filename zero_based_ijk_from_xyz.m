function result = zero_based_ijk_from_xyz(xyz, origin, spacing)
    % xyz is n x 3, in um
    % offset is 1 x 3, in um
    % spacing is 1 x 3, in um
    xyz_centered = bsxfun(@minus, xyz, origin) ;
    ijk_zero_based_continuous = bsxfun(@rdivide, xyz_centered, spacing) ;
    result = floor(ijk_zero_based_continuous) ;
end
