function result = ijk1_from_jaws_xyz(xyz, origin, spacing)
    % xyz is n x 3, in um
    % offset is 1 x 3, in um
    % spacing is 1 x 3, in um
    ijk0 = ijk0_from_jaws_xyz(xyz, origin, spacing) ;
    result = ijk0 + 1 ;
end
