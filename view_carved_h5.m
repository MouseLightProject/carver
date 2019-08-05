octree_path = '/nrs/mouselight/SAMPLES/2018-10-01' ;
h5_file_path = '/home/taylora/cache/gt/2018-10-01/navigator-output/test-swcs-carved.h5' ;
%offset = [75445.691517998261       16536.315191998805       36051.105973001562] ;  % soma 1
offset = [75885.516008994338       15382.253859001321        35947.93727700039] ;  % soma 2
%sample_point_xyz_centered = [2074.0422090017528      -684.41304799880709    1333.5478669984368] ;  % soma 1 
sample_point_xyz_centered = [1157.5993180056539        11.56324699867946       1683.9443079996054] ;  % soma 2
%sample_point_xyz_centered = [738.19657500174071      -1424.7829449988058       702.41362699843739] ;
sample_point_xyz = sample_point_xyz_centered + offset

[origin, spacing] = load_transform_txt(fullfile(octree_path, 'transform.txt')) ;
sample_point_ijk = one_based_ijk_from_xyz(sample_point_xyz, origin, spacing)
%half_diagonal = [512 512 128] ;
half_diagonal = 1024 * [1 1 1/4] ;
stack_lower_corner = sample_point_ijk - half_diagonal
stack_upper_corner = sample_point_ijk + (half_diagonal-1)
stack_shape = 2*half_diagonal ;

info  = h5info(h5_file_path, '/volume')
raw_stack_4d = h5read(h5_file_path, '/volume', fliplr([stack_lower_corner 1]), fliplr([stack_shape 1])) ;  % stupid color channels
stack_4d = permute(raw_stack_4d, fliplr(1:4)) ;
shape_4d = size(stack_4d) ;
shape = shape_4d(1:3) ;
stack = reshape(stack_4d, shape) ;

volumeViewer(stack) ;
