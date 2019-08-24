function dump_string_to_file(file_name, str)
    fid = fopen(file_name, 'wt') ;
    fprintf(fid, '%s', str) ;
    fclose(fid) ;    
end