function [G, subs] = load_skeleton_graph_from_txt_files(skeleton_folder_name, sample_stack_shape)
    % sample_stack_shape is in xyz order.  This is the shape of the
    % rendered sample, at the highest resolution in the 'octree'.

    % Verify the folder exists
    if ~exist(skeleton_folder_name, 'file') ,
        error('The skeletonization output folder %s does not exist', skeleton_folder_name) ;
    end
    
    folder_listing = dir(fullfile(skeleton_folder_name,'*.txt'));
    % check the format for a non zero file
    txt_file_names = {folder_listing.name} ;
    txt_file_sizes = [folder_listing.bytes];
    
    % Figure out the number of columns in each txt file
    index_of_first_nonempty_file = find(txt_file_sizes>0,1) ;
    fid = fopen(fullfile(skeleton_folder_name, txt_file_names{index_of_first_nonempty_file}));
    tline = fgetl(fid);
    fclose(fid);
    tlines = strsplit(tline,' ');
    column_count = size(tlines,2); % first 2 are indicies, rest are weights
    
    % Determine the format based on the column count
    format = ['%f %f', sprintf('%s',repmat(' %f',ones(1,column_count-2)))];

    txt_file_count = length(txt_file_names) ;
    pairs = cell(1,txt_file_count) ;
    use_all_cores() ;
    parfor idx = 1:txt_file_count
        if txt_file_sizes(idx) == 0 ,
            continue
        end
        % read text file line by line
        fid = fopen(fullfile(skeleton_folder_name, txt_file_names{idx})) ;
        tmp = textscan(fid, format) ;
        fclose(fid);
        pairs{idx} = cat(2,tmp{:})';
    end
    raw_edges = [pairs{:}]' ; 

%     [keepthese,ia,ic] = unique(edges(:, [1 2])) ;
    [keepthese, ~, ic] = unique(raw_edges(:, [1 2])) ;
    [subs(:,1),subs(:,2),subs(:,3)] = ind2sub(sample_stack_shape([1 2 3]), keepthese) ;
    edges = reshape(ic, [], 2) ;
%     weights = edges(ia,3:end);
    if isempty(edges) ,
        return
    end
%     if isempty(weights) ,
%         weights = ones(size(edges_,1),1);
%     end

    selfinds = find((edges(:,1)==edges(:,2))) ;
    if ~isempty(selfinds) ,
        edges(selfinds,:) = [] ;
    end
    
    A_raw = sparse(edges(:,1),edges(:,2),1,max(edges(:)),max(edges(:)));
    A = max(A_raw',A_raw);
    G = graph(A) ;
end
