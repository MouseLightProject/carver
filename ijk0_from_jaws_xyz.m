function result = ijk0_from_jaws_xyz(xyz, origin, spacing)
    % xyz is n x 3, in um
    % offset is 1 x 3, in um
    % spacing is 1 x 3, in um
    xyz_centered = xyz -origin ;
    ijk0_continuous = xyz_centered ./ spacing ;
    result = floor(ijk0_continuous) ;
end
