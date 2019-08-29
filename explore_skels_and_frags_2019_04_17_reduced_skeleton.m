sample_date = '2019-04-17' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
%output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/%s/carver', sample_date) ;
output_folder_path = fullfile(this_folder_path, sprintf('%s-reduced-skeleton', sample_date)) ;
%trees_folder_path = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/build-brain-output/full', sample_date) ;  % these are .swcs for this sample
trees_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/scripts/build-brain/%s-reduced-skeleton-output', sample_date) ;
sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
skeleton_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/scripts/build-brain/%s-reduced-skeleton', sample_date) ;



% Get the metadata for the rendered sample
[shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(sample_folder_path);

% % Load the fragments from disk                                                
% all_fragment_xyzs_in_erhan_coords = ...
%     compute_or_read_from_memo(output_folder_path, ...
%                               'all_fragment_xyzs_in_erhan_coords',  ...
%                               @()(read_all_fragment_centerpoints(fragments_folder_path))) ;
% 
% % Check that the fragment nodes are on the grid as we expect
% %     all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
% %     all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
% %     fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
% %     assert( all( all( abs(fragment_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
% fraction_of_fragment_nodes_at_voxel_centers = fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(all_fragment_xyzs_in_erhan_coords, spacing, origin)
% assert( fraction_of_fragment_nodes_at_voxel_centers == 1) ;
% 
% % Convert fragment coords to JaWS coords
% all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
% %all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 1.5 ;  % seems like maybe this is correct?
% all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
% all_fragment_xyzs = jaws_origin + spacing .* (all_fragment_ijks-1) ;  % um, n x 3

% Get the skeleton graph, either from .txt files or from mat          
if ~exist(output_folder_path, 'file') ,
    mkdir(output_folder_path) ;
end
graph_file_path = fullfile(output_folder_path, 'skeleton-as-graph.mat') ;
if exist(graph_file_path, 'file') ,
    load(graph_file_path, 'skeleton_graph', 'skeleton_ijks') ;
else        
    [skeleton_graph, skeleton_ijks] = load_skeleton_graph_from_txt_files(skeleton_folder_path, shape_xyz) ;  % these ijks are also one-based
    save(graph_file_path, 'skeleton_graph','skeleton_ijks', '-v7.3') ;
end

% % What fraction of fragment points are also skeleton points, looking at ijk
% % coords
% is_fragment_a_skeleton_point_from_ijk = is_point_drawn_from_pool(all_fragment_ijks, skeleton_ijks, [1 1 1]) ;
% fraction_fragment_points_that_are_skeleton_points_from_ijk = mean(is_fragment_a_skeleton_point_from_ijk)
% 
% % What if you shift the fragment points by 1 in each dim?
% is_shifted_fragment_a_skeleton_point_from_ijk = is_point_drawn_from_pool(all_fragment_ijks+1, skeleton_ijks, [1 1 1]) ;
% fraction_shifted_fragment_points_that_are_skel_points_from_ijk = mean(is_shifted_fragment_a_skeleton_point_from_ijk)


skeleton_xyzs = jaws_origin + spacing .* (skeleton_ijks-1) ;  % um, n x 3
    % skeleton_ijks are 1-based, we want to map them to voxel centers
    % in JaWS coordinates

% % Write the skeleton to json file, for Will & Jan
% graph_as_json_file_path = fullfile(output_folder_path, 'skeleton-as-graph.json') ;
% if ~exist(graph_as_json_file_path, 'file') ,
%     json_ready_skeleton = struct() ;
%     json_ready_skeleton.ijks = skeleton_ijks-1 ;  % convert to zero-based indexing
%     edges_with_third_col = table2array(skeleton_graph.Edges)-1 ;  % convert to zero-based indexing        
%     json_ready_skeleton.edges = edges_with_third_col(:,1:2) ;
%     json_ready_skeleton.jaws_origin = jaws_origin ;
%     json_ready_skeleton.spacing = spacing ;
%     json_ready_skeleton.origin = origin ;
%     json_ready_skeleton.shape_xyz = shape_xyz ;        
%     json_ready_skeleton.xyz_from_ijk_formula = 'xyzs = jaws_origin + spacing .* ijks' ;        
%     json_text = jsonencode(json_ready_skeleton) ;
%     dump_string_to_file(graph_as_json_file_path, json_text) ;
% end    

% % What fraction of fragment points are also skeleton points
% is_fragment_a_skeleton_point = is_point_drawn_from_pool(all_fragment_xyzs, skeleton_xyzs, spacing) ;
% fraction_fragment_points_that_are_skeleton_points = mean(is_fragment_a_skeleton_point)

% Make fragment points, skeleton point indices 
%fragments_kd_tree  = KDTreeSearcher(all_fragment_xyzs) ;        
skeleton_kd_tree = KDTreeSearcher(skeleton_xyzs) ;        

% %%
% % Pick a random fragment point, plot all the fragment and skeleton points
% % nearby.
% fragment_node_count = size(all_fragment_xyzs, 1)
% fragment_node_index = randi(fragment_node_count, 1)
% test_xyz = all_fragment_xyzs(fragment_node_index, :)
% 
% r = 10 ;  % um
% indices_of_nearby_fragment_points_as_cell = fragments_kd_tree.rangesearch(test_xyz, r) ;
% indices_of_nearby_fragment_points = indices_of_nearby_fragment_points_as_cell{1} ;
% xyzs_of_nearby_fragment_points = all_fragment_xyzs(indices_of_nearby_fragment_points, :) ;
% indices_of_nearby_skeleton_points_as_cell = skeleton_kd_tree.rangesearch(test_xyz, r) ;
% indices_of_nearby_skeleton_points = indices_of_nearby_skeleton_points_as_cell{1} ;
% xyzs_of_nearby_skeleton_points = skeleton_xyzs(indices_of_nearby_skeleton_points, :) ;
% 
% f = figure('Color', 'w') ;
% a = axes(f) ;
% l1 = line_in_3d_bang(a, xyzs_of_nearby_fragment_points-test_xyz, 'LineStyle', 'none', 'Marker', '.', 'Color', 'b') ;
% l2 = line_in_3d_bang(a, xyzs_of_nearby_skeleton_points-test_xyz, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r') ;
% l0 = line_in_3d_bang(a, test_xyz-test_xyz, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 12, 'Color', 'b') ;
% a.DataAspectRatio=[1 1 1] ;
% view(3) ;
% axis vis3d
% xlabel('x (um)') ;
% ylabel('y (um)') ;
% zlabel('z (um)') ;
% drawnow;
% xlim(1.01*r*[-1 +1])
% ylim(1.01*r*[-1 +1])
% zlim(1.01*r*[-1 +1])



% %%
% desired_sample_count = 100 ;
% [sample_stacks, sample_radius] = sample_rendered_data(sample_folder_path, all_fragment_ijks, desired_sample_count) ;
% mean_sample_stack = mean(sample_stacks, 4) ;
% mean_sample_stack_xy_slice = mean_sample_stack(:,:,sample_radius(3)+1) ;
% f = figure('Color', 'w') ;
% a = axes(f, 'DataAspectRatio', [1 1 1]) ;
% imagesc(a, sample_radius(1)*[-1 +1], sample_radius(1)*[-1 +1], mean_sample_stack_xy_slice) ;
% axis equal
% axis tight
% line(a, 'XData', 0, 'YData', 0, 'Marker', '.', 'MarkerSize', 12) ;
% title('Fragments') ;

% %%
% desired_sample_count = 200 ;
% [sample_stacks, sample_radius] = sample_rendered_data(sample_folder_path, skeleton_ijks, desired_sample_count) ;
% mean_sample_stack = mean(sample_stacks, 4) ;
% mean_sample_stack_xy_slice = mean_sample_stack(:,:,sample_radius(3)+1) ;
% f = figure('Color', 'w') ;
% a = axes(f, 'DataAspectRatio', [1 1 1]) ;
% imagesc(a, sample_radius(1)*[-1 +1], sample_radius(1)*[-1 +1], mean_sample_stack_xy_slice) ;
% axis equal
% axis tight
% line(a, 'XData', 0, 'YData', 0, 'Marker', '.', 'MarkerSize', 12) ;
% title('Skeleton') ;

% sample_stack_xy_slices = squeeze(sample_stacks(:,:,sample_radius(3)+1,:)) ;
% f = figure('Color', 'w') ;
% a = axes(f, 'DataAspectRatio', [1 1 1]) ;
% imagesc(a, sample_radius(1)*[-1 +1], sample_radius(1)*[-1 +1], sample_stack_xy_slices(:,:,60)) ;
% axis equal
% axis tight
% line(a, 'XData', 0, 'YData', 0, 'Marker', '.', 'MarkerSize', 12) ;
% title('Skeleton') ;



%%
% Load the tree points from disk                                                
all_tree_xyzs_in_erhan_coords = ...
    compute_or_read_from_memo(output_folder_path, ...
                              'all_tree_xyzs_in_erhan_coords',  ...
                              @()(read_all_centerpoints_in_folder(trees_folder_path, '.mat'))) ;

% Check that the tree nodes are on the grid as we expect
%     all_tree_ijks_unrounded = (all_tree_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
%     all_tree_ijks = round(all_tree_ijks_unrounded) ;  % these should be one-based ijks
%     tree_offsets = all_tree_ijks_unrounded - all_tree_ijks ;
%     assert( all( all( abs(tree_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
fraction_of_tree_nodes_at_voxel_centers = fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(all_tree_xyzs_in_erhan_coords, spacing, origin)

% Convert tree coords to JaWS coords
all_tree_ijks_unrounded = (all_tree_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
all_tree_ijks = round(all_tree_ijks_unrounded) ;  % these should be one-based ijks
all_tree_xyzs = jaws_origin + spacing .* (all_tree_ijks-1) ;  % um, n x 3

% What fraction of tree points are also skeleton points, looking at ijk
% coords
is_tree_a_skeleton_point_from_ijk = is_point_drawn_from_pool(all_tree_ijks, skeleton_ijks, [1 1 1]) ;
fraction_tree_points_that_are_skeleton_points_from_ijk = mean(is_tree_a_skeleton_point_from_ijk)

% What if you shift the tree points by 1 in each dim?
is_shifted_tree_a_skeleton_point_from_ijk = is_point_drawn_from_pool(all_tree_ijks+1, skeleton_ijks, [1 1 1]) ;
fraction_shifted_tree_points_that_are_skel_points_from_ijk = mean(is_shifted_tree_a_skeleton_point_from_ijk)



