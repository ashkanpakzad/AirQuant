# TODO: Very sloppy imports - need to fix
from funcs_br import *
import sknw.sknw.sknw as sknw # for some reason the normal import does not work

dat_path = 'C:/PhD/test_data'

t = time.time()

# Read with itk.
image = itk.imread(dat_path + '/case21.nii')
seg = itk.imread(dat_path + '/case21_seg.nii.gz')

# Run skeletonisation filter - good performance ~ 10s
# Based on Homann H., Implementation of a 3D thinning algorithm, The Insight Journal - 2007
# ensures 26 connectivity
print('Performing ITK skeletonization...')
asf = itk.BinaryThinningImageFilter3D.New(seg)

# ~ 3D distance transform - fits a sphere - can use as guideline
print('Performing ITK 3D distance transform...')
thickness_map = itk.MedialThicknessImageFilter3D.New(seg)

# Saving the itk file
# itk.imwrite(thickness_map, dat_path + '/results/itk_thickmap.nii.gz')

# Copy of itk.Image to np array, data is copied - takes ~2.5 min
print('Converting skeleton to np...')
skeleton = itk.array_from_image(asf)

# Build the graph - this actually works but raises warnings because of jit
# TODO: edit sknw to fix JIT speed-up | update or remove JIT | suppress warnings
print('Converting to graph with sknw...')
graph = sknw.build_sknw(skeleton)

# 3D Spline testing - for 1 edge only
edge56 = graph.edges[(5,6)]['pts']
print('Fitting splines...')

# Use 5 point smoothing (s=5) as Kin did; default is s = m - sqrt(2*m) with m = # datapoints
# TODO: check whether smoothing in splprep has same behaviour as in Matlab
tck, u = interpolate.splprep([edge56[...,0], edge56[...,1], edge56[...,2]], s=5)
x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
u_fine = np.linspace(0,1,100)
x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

elapsed = time.time() - t
print('Time elapsed = ', elapsed)


# Plotting
# Node coordinates for plotting
graph_vertices = np.empty(3,)
for i in range(len(graph.nodes)):
    # this produces 'average' coords, use 'pts' to see when multiple pts belong to a node - e.g. pts around bifurcation
    chunk = graph.nodes[i]['o']
    graph_vertices = np.vstack([graph_vertices, chunk])

plot_3d(skeleton, graph_vertices,0)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(edge56[...,0],edge56[...,1], edge56[...,2], s=5, c='red', marker='o')
ax.scatter(x_fine, y_fine, z_fine, c='blue')


### Unused code below this ###

# segmentation = nib.load(os.path.join(dat_path, 'case21_seg.nii'))
# seg = segmentation.get_fdata()
# raw_image = nib.load(os.path.join(dat_path, 'case21.nii'))
# img = raw_image.get_fdata()

# View only of itk.Image, data is not copied
# np_view = itk.array_view_from_image(asf)

# Write with ITK - preserves orientation
# itk.imwrite(asf, dat_path + '/results/itk_skel.nii.gz')

# Attempt saving with nib - saves 90 degrees rotated, need to fix
# segmentation = nib.load(os.path.join(dat_path, 'case21_seg.nii'))
# savepath = os.path.join(dat_path, 'results')
# save_image(splined, segmentation.affine, segmentation.header, savepath+'/splined_post_itk')
