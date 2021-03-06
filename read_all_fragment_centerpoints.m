function result = read_all_fragment_centerpoints(fragments_folder_path)
    if ~exist(fragments_folder_path, 'file') ,
        error('The folder %s does not exist', fragments_folder_path) ;
    end       
    fragment_swc_template = fullfile(fragments_folder_path, '*.swc') ;
    fragment_file_names = simple_dir(fragment_swc_template) ;
    fragment_count = length(fragment_file_names)
    raw_result = zeros(10e6, 3) ;
    centerpoint_count = 0 ;
    pbo = progress_bar_object(fragment_count) ;
    for fragment_index = 1 : fragment_count ,
        fragment_file_name = fragment_file_names{fragment_index} ;
        fragment_file_path = fullfile(fragments_folder_path, fragment_file_name) ;
        swc_data = load_swc(fragment_file_path) ;
        fragment_centerpoints = swc_data(:,3:5) ;  % each row an xyz triple, in um
        if any(min(fragment_centerpoints, 1)<0) ,
            error('badness!') ;
        end        
        fragment_centerpoint_count = size(fragment_centerpoints, 1) ;
        raw_result(centerpoint_count+1:centerpoint_count+fragment_centerpoint_count, :) = fragment_centerpoints ;
        centerpoint_count = centerpoint_count + fragment_centerpoint_count ;        
        pbo.update(fragment_index) ;
    end    
    result = raw_result(1:centerpoint_count,:) ;
end
