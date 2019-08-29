function result = read_all_centerpoints_in_folder(folder_name, extension)
    if ~exist(folder_name, 'file') ,
        error('The folder %s does not exist', folder_name) ;
    end       
    file_name_template = fullfile(folder_name, sprintf('auto*%s', extension)) ;
    file_names = simple_dir(file_name_template) ;
    file_count = length(file_names)
    raw_result = zeros(10e6, 3) ;
    total_centerpoint_count = 0 ;
    pbo = progress_bar_object(file_count) ;
    for file_index = 1 : file_count ,
        file_name = file_names{file_index} ;
        file_path = fullfile(folder_name, file_name) ;
        if isequal(extension, '.swc') ,
            swc_as_array = load_swc(file_path) ;
            centerpoints = swc_as_array(:,3:5) ;  % each row an xyz triple, in um
        else
            s = load(file_path, '-mat') ;
            swc_as_struct = s.outtree ;
            centerpoints = [swc_as_struct.X swc_as_struct.Y swc_as_struct.Z] ;
        end
        if any(min(centerpoints, 1)<0) ,
            error('badness!') ;
        end        
        file_centerpoint_count = size(centerpoints, 1) ;
        raw_result(total_centerpoint_count+1:total_centerpoint_count+file_centerpoint_count, :) = centerpoints ;
        total_centerpoint_count = total_centerpoint_count + file_centerpoint_count ;        
        pbo.update(file_index) ;
    end    
    result = raw_result(1:total_centerpoint_count,:) ;
end
