function save_swc(file_name, swc_data, neuron_name)
    was_neuron_name_given =  ( exist('neuron_name', 'var') && ~isempty(neuron_name) ) ;
    
    xyz = swc_data(:,3:5) ;
    offset = mean(xyz,1) ;
    xyz_centered = bsxfun(@minus, xyz, offset) ;
    centered_swc_data = swc_data ;
    centered_swc_data(:,3:5) = xyz_centered ;
    
    fid = fopen(file_name,'wt') ;
    fprintf(fid, '# Generated by save_swc.m\n') ;
    fprintf(fid, '# OFFSET %24.17g %24.17g %24.17g\n', offset) ;
    fprintf(fid, '# COLOR 1.000000,1.000000,1.000000\n', offset) ;
    if was_neuron_name_given ,
        fprintf(fid, '# NAME %s.swc\n', neuron_name) ;
    end
    fprintf(fid,'%d %d %24.17g %24.17g %24.17g %d %d\n',centered_swc_data') ;
    fclose(fid) ;
end
