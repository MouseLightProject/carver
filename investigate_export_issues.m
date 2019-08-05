original_neuron = load_swc('test-tree/auto111_cc-0231.swc')
exported_neuron = load_swc('auto111_cc-0231-as-exported.swc')
exported_neuron - original_neuron  % all off by same amount
mean(original_neuron)
diff(original_neuron,1)  
diff(exported_neuron,1)  % same (up to floating-point error) as last value
[origin, spacing] = load_transform_txt('/nrs/mouselight/SAMPLES/2018-10-01/transform.txt')
origin = [69445.01944245792, 12917.29486937314, 30198.96941474185] ;  % origin from Workstation
spacing(3)/spacing(1)  % about 4
spacing(3)/spacing(2)  % about 4
original_neuron_xyz = original_neuron(:,3:5)
exported_neuron_xyz = exported_neuron(:,3:5)
diff(original_neuron_xyz)
dr_original = diff(original_neuron_xyz)
bsxfun(@times, original_neuron_xyz-origin, 1./spacing)  % half-integers
bsxfun(@times, exported_neuron_xyz-origin, 1./spacing)  % not half-integers
exported_centroid = mean(exported_neuron_xyz)
bsxfun(@times, exported_centroid-origin, 1./spacing)  % not half-integers, or integers...
exported_neuron_xyz - original_neuron_xyz  % all columns the same
original_centroid = mean(original_neuron_xyz)
exported_centroid-original_centroid
