close all
clear

sample_date = '2018-07-02' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
%output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/%s/carver', sample_date) ;
output_folder_path = fullfile(this_folder_path, sprintf('%s-v03', sample_date)) ;
fragments_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/prob0_swcs/frags', sample_date) ;
consensus_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;
sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
skeleton_folder_path =  sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;
xyz_of_interest = [78607.5, 19303.7, 31703.7] ;  % um, JAWS coords


% Get the metadata for the rendered sample
[shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(sample_folder_path);

% Load the consensus neuron
consensus_swc_path = '/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/2018-07-02/G-002/consensus/2018-07-02_G-002_consensus.swc' ;
consensus_neuron_as_swc_array = load_swc(consensus_swc_path) ;
consensus_xyzs = consensus_neuron_as_swc_array(:,3:5) ;

% Extract a stretch of consensus
[min_dist, i_min_dist] = min(sqrt(sum((consensus_xyzs-xyz_of_interest).^2, 2)))
nearby_consensus_xyzs = consensus_xyzs(i_min_dist-5:i_min_dist+5,:) ;

% Load the refined swc
refined_swc_path = '/groups/mousebrainmicro/mousebrainmicro/scripts/carver/2018-07-02-v03/augmented-with-skeleton-nodes-as-swcs/G-002.swc' ;
refined_neuron_as_swc_array = load_swc(refined_swc_path) ;
refined_xyzs = refined_neuron_as_swc_array(:,3:5) ;

% Extract a stretch of refined
[min_dist, i_min_dist] = min(sqrt(sum((refined_xyzs-xyz_of_interest).^2, 2)))
nearby_refined_xyzs = refined_xyzs(i_min_dist-100:i_min_dist+100,:) ;



ijk1_of_interest = round((xyz_of_interest - jaws_origin) ./ spacing) + 1 ;
sample_shape = [64 80 6] ;
sample_offset_ijk1 = ijk1_of_interest - round(sample_shape/2) ;
sample_stack = get_mouselight_rendered_substack(sample_folder_path, 0, sample_offset_ijk1, sample_shape) ;
sample_stack_mip = max(sample_stack, [], 3) ;

sample_high_corner_ijk1 = sample_offset_ijk1 + sample_shape - 1 ;
sample_offset_xyz = jaws_xyz_from_ijk1(sample_offset_ijk1, jaws_origin, spacing) ;
sample_high_corner_xyz = jaws_xyz_from_ijk1(sample_high_corner_ijk1, jaws_origin, spacing) ;

f = figure('Color', 'w') ;
a = axes(f, 'DataAspectRatio', [1 1 1]) ;
imagesc(a, ...
        [sample_offset_xyz(1) sample_high_corner_xyz(1)], ...
        [sample_offset_xyz(2) sample_high_corner_xyz(2)], ...
        sample_stack_mip) ;
f.Colormap = green(256) ;    
a.DataAspectRatio = [1 1 1] ;
xlim([sample_offset_xyz(1)-spacing(1)/2 sample_high_corner_xyz(1)+spacing(1)/2]) ;
ylim([sample_offset_xyz(2)-spacing(2)/2 sample_high_corner_xyz(2)+spacing(2)/2]) ;
a.XTick = [] ;
a.YTick = [] ;

consensus_line = line_in_3d_bang(a, nearby_consensus_xyzs, 'Marker', '.', 'MarkerSize', 3*6, 'Color', 'r') ;
refined_line = line_in_3d_bang(a, nearby_refined_xyzs, 'Marker', '.', 'MarkerSize', 3*6, 'Color', [0 0.5 1]) ;

set_figure_size([8 9]) ;
set_figure_to_wysiwyg_printing(f) ;
f.Name = 'example-refinement' ;
print_pdf(f) ;


% Check the alignment of refined neuron xyz points with the imagery
%%
refined_ijk1s = ijk1_from_jaws_xyz(refined_xyzs, jaws_origin, spacing) ;
desired_sample_count = 100 ;
[sample_stacks, sample_radius] = sample_rendered_data(sample_folder_path, refined_ijk1s, desired_sample_count) ;
mean_sample_stack = mean(sample_stacks, 4) ;
mean_sample_stack_xy_slice = mean_sample_stack(:,:,sample_radius(3)+1) ;
f = figure('Color', 'w') ;
a = axes(f, 'DataAspectRatio', [1 1 1]) ;
imagesc(a, sample_radius(1)*[-1 +1], sample_radius(1)*[-1 +1], mean_sample_stack_xy_slice) ;
axis equal
axis tight
line(a, 'XData', 0, 'YData', 0, 'Marker', '.', 'MarkerSize', 12) ;
title('Refined points') ;







% Load the fragments from disk                                                
all_fragment_xyzs_in_erhan_coords = ...
    compute_or_read_from_memo(output_folder_path, ...
                              'all_fragment_xyzs_in_erhan_coords',  ...
                              @()(read_all_fragment_centerpoints(fragments_folder_path))) ;

% Check that the fragment nodes are on the grid as we expect
%     all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
%     all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
%     fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
%     assert( all( all( abs(fragment_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
fraction_of_fragment_nodes_at_voxel_centers = fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(all_fragment_xyzs_in_erhan_coords, spacing, origin)
assert( fraction_of_fragment_nodes_at_voxel_centers == 1) ;

% Convert fragment coords to JaWS coords
all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
%all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 1.5 ;  % seems like maybe this is correct?
all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
all_fragment_xyzs = jaws_origin + spacing .* (all_fragment_ijks-1) ;  % um, n x 3

% Get the skeleton graph, either from .txt files or from mat                              
graph_file_path = fullfile(output_folder_path, 'skeleton-as-graph.mat') ;
if exist(graph_file_path, 'file') ,
    load(graph_file_path, 'skeleton_graph', 'skeleton_ijks') ;
else        
    [skeleton_graph, skeleton_ijks] = load_skeleton_graph_from_txt_files(skeleton_folder_path, shape_xyz) ;  % these ijks are also one-based
    save(graph_file_path, 'skeleton_graph','skeleton_ijks', '-v7.3') ;
end

% What fraction of fragment points are also skeleton points, looking at ijk
% coords
is_fragment_a_skeleton_point_from_ijk = is_point_drawn_from_pool(all_fragment_ijks, skeleton_ijks, [1 1 1]) ;
fraction_fragment_points_that_are_skeleton_points_from_ijk = mean(is_fragment_a_skeleton_point_from_ijk)

% What if you shift the fragment points by 1 in each dim?
is_shifted_fragment_a_skeleton_point_from_ijk = is_point_drawn_from_pool(all_fragment_ijks+1, skeleton_ijks, [1 1 1]) ;
fraction_shifted_fragment_points_that_are_skeleton_points_from_ijk = mean(is_shifted_fragment_a_skeleton_point_from_ijk)


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

% What fraction of fragment points are also skeleton points
is_fragment_a_skeleton_point = is_point_drawn_from_pool(all_fragment_xyzs, skeleton_xyzs, spacing) ;
fraction_fragment_points_that_are_skeleton_points = mean(is_fragment_a_skeleton_point)

% Make fragment points, skeleton point indices 
fragments_kd_tree  = KDTreeSearcher(all_fragment_xyzs) ;        
skeleton_kd_tree = KDTreeSearcher(skeleton_xyzs) ;        

%%
% Pick a random fragment point, plot all the fragment and skeleton points
% nearby.
fragment_node_count = size(all_fragment_xyzs, 1)
fragment_node_index = randi(fragment_node_count, 1)
test_xyz = all_fragment_xyzs(fragment_node_index, :)

r = 10 ;  % um
indices_of_nearby_fragment_points_as_cell = fragments_kd_tree.rangesearch(test_xyz, r) ;
indices_of_nearby_fragment_points = indices_of_nearby_fragment_points_as_cell{1} ;
xyzs_of_nearby_fragment_points = all_fragment_xyzs(indices_of_nearby_fragment_points, :) ;
indices_of_nearby_skeleton_points_as_cell = skeleton_kd_tree.rangesearch(test_xyz, r) ;
indices_of_nearby_skeleton_points = indices_of_nearby_skeleton_points_as_cell{1} ;
xyzs_of_nearby_skeleton_points = skeleton_xyzs(indices_of_nearby_skeleton_points, :) ;

f = figure('Color', 'w') ;
a = axes(f) ;
l1 = line_in_3d_bang(a, xyzs_of_nearby_fragment_points-test_xyz, 'LineStyle', 'none', 'Marker', '.', 'Color', 'b') ;
l2 = line_in_3d_bang(a, xyzs_of_nearby_skeleton_points-test_xyz, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r') ;
l0 = line_in_3d_bang(a, test_xyz-test_xyz, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 12, 'Color', 'b') ;
a.DataAspectRatio=[1 1 1] ;
view(3) ;
axis vis3d
xlabel('x (um)') ;
ylabel('y (um)') ;
zlabel('z (um)') ;
drawnow;
xlim(1.01*r*[-1 +1])
ylim(1.01*r*[-1 +1])
zlim(1.01*r*[-1 +1])



%%
desired_sample_count = 100 ;
[sample_stacks, sample_radius] = sample_rendered_data(sample_folder_path, all_fragment_ijks, desired_sample_count) ;
mean_sample_stack = mean(sample_stacks, 4) ;
mean_sample_stack_xy_slice = mean_sample_stack(:,:,sample_radius(3)+1) ;
f = figure('Color', 'w') ;
a = axes(f, 'DataAspectRatio', [1 1 1]) ;
imagesc(a, sample_radius(1)*[-1 +1], sample_radius(1)*[-1 +1], mean_sample_stack_xy_slice) ;
axis equal
axis tight
line(a, 'XData', 0, 'YData', 0, 'Marker', '.', 'MarkerSize', 12) ;
title('Fragments') ;

%%
desired_sample_count = 100 ;
[sample_stacks, sample_radius] = sample_rendered_data(sample_folder_path, skeleton_ijks, desired_sample_count) ;
mean_sample_stack = mean(sample_stacks, 4) ;
mean_sample_stack_xy_slice = mean_sample_stack(:,:,sample_radius(3)+1) ;
f = figure('Color', 'w') ;
a = axes(f, 'DataAspectRatio', [1 1 1]) ;
imagesc(a, sample_radius(1)*[-1 +1], sample_radius(1)*[-1 +1], mean_sample_stack_xy_slice) ;
axis equal
axis tight
line(a, 'XData', 0, 'YData', 0, 'Marker', '.', 'MarkerSize', 12) ;
title('Skeleton') ;







