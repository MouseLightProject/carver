function makecellfun(lambda, cell_array, file_name_template)
    n = length(cell_array) ;
    for i = 1:n ,
        file_name = sprintf(file_name_template, i) ;
        if ~exist(file_name, 'file') ,
            feval(lambda, cell_array{i}, file_name) ;
        end
    end
end
