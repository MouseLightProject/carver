function use_all_cores()
    poolobj = gcp('nocreate');  % Get the current pool, but don't create one if it doesn't exist yet
    if isempty(poolobj) ,
        % If no pool exists, create one with as many processes as there are
        % physical cores.
        parpool(feature('numcores'))
    else
        % If a pool already exists, use that one.
    end
end
