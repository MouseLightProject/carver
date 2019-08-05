function result = one_based_ijk_from_xyz(xyz, origin, spacing)
    % xyz is n x 3, in um
    % offset is 1 x 3, in um
    % spacing is 1 x 3, in um
    ijk_zero_based = zero_based_ijk_from_xyz(xyz, origin, spacing) ;
    result = ijk_zero_based + 1 ;
end
