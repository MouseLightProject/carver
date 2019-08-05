function result = read_all_fragment_centerpoints(fragments_folder_path)
    fragment_swc_template = fullfile(fragments_folder_path, '*.swc') ;
    fragment_file_names = simple_dir(fragment_swc_template) ;
    fragment_count = length(fragment_file_names)
    raw_result = zeros(10e6, 3) ;
    centerpoint_count = 0 ;
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
        if mod(fragment_index, 100) == 0 ,
            fprintf('.') ;
            if mod(fragment_index, 10000) == 0 ,
                fprintf('\n') ;
            end
        end
    end
    % Print a newline, unless we just printed one
    if mod(fragment_count, 100000) ~= 0 ,
        fprintf('\n') ;
    end
    
    result = raw_result(1:centerpoint_count,:) ;
end
