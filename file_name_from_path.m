function result = file_name_from_path(path) 
    [~, base_name, ext] = fileparts(path) ;
    result = horzcat(base_name, ext) ;
end
