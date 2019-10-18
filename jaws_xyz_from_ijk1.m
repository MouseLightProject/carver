function result = jaws_xyz_from_ijk1(ijk1, origin, spacing)
    % xyz is n x 3, in um
    % offset is 1 x 3, in um
    % spacing is 1 x 3, in um
    ijk0 = ijk1 - 1 ;
    result = origin + spacing .* ijk0 ;
end
