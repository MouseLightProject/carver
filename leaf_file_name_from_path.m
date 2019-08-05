function result = leaf_file_name_from_path(path) 
    [~,base,ext] = fileparts(path) ;
    result = horzcat(base, ext) ;
end
