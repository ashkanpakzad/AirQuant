# TODO: Very sloppy imports - need to fix
from funcs_br import *
# import sknw.sknw.sknw as sknw # for some reason the normal import does not work
import sknw
import scipy

# TODO: crossplatform dataloading
dat_path = '/Users/apakz/PhD/bojidar_data/case21'

t = time.time()

### Read with itk.
image = itk.imread(dat_path + '/case21.nii')
seg = itk.imread(dat_path + '/case21_seg.nii.gz')

### Run skeletonisation filter - good performance ~ 10s
# Based on Homann H., Implementation of a 3D thinning algorithm, The Insight Journal - 2007
# ensures 26 connectivity
# TODO: consider providing utility for lesser points of connectivity
print('Performing ITK skeletonization...')
asf = itk.BinaryThinningImageFilter3D.New(seg)

### ~ 3D distance transform - fits a sphere - can use as guideline
print('Performing ITK 3D distance transform...')
thickness_map = itk.MedialThicknessImageFilter3D.New(seg)

# Saving the itk file
# itk.imwrite(thickness_map, dat_path + '/results/itk_thickmap.nii.gz')
# Copy of itk.Image to np array, data is copied - takes ~2.5 min -- took my lil laptop 9 mins :(
print('Converting skeleton to np...')
skeleton = itk.array_from_image(asf)

# Build the graph - this actually works but raises warnings because of jit
# TODO: edit sknw to fix JIT speed-up | update or remove JIT | suppress warnings
print('Converting to graph with sknw...')
graph = sknw.build_sknw(skeleton)

### Plot skeleton and graph
# Node coordinates for plotting
graph_vertices = np.empty(3, )
for i in range(len(graph.nodes)):
    # this produces 'average' coords, use 'pts' to see when multiple pts belong to a node - e.g. pts around bifurcation
    chunk = graph.nodes[i]['o']
    graph_vertices = np.vstack([graph_vertices, chunk])

plot_3d(skeleton, graph_vertices, 0)

### 3D Spline testing - for 1 edge only
edge56 = graph.edges[(5, 6)]['pts']
print('Fitting splines...')

# Use 5 point smoothing (s=5) as Kin did; default is s = m - sqrt(2*m) with m = # datapoints
# TODO: check whether smoothing in splprep has same behaviour as in Matlab
tck, u = interpolate.splprep([edge56[..., 0], edge56[..., 1], edge56[..., 2]], s=5)
# tck = spline object, u = parametrisation along spline.
# x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)

### Plot spline along data
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(edge56[..., 0], edge56[..., 1], edge56[..., 2], s=5, c='red', marker='o')
ax.scatter(x_fine, y_fine, z_fine, c='blue')

## Compute tangents along spline
# TODO: consider how to construct spline interpolation interval.
u_fine = np.linspace(0, 1, 17)  # interpolate at each point in edge for now
x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)  # get spline coords
# TODO: work out how to get the tangent!
def central_diff(vector, point):
    return (vector[point + 1] - vector[point - 1]) / 2

spline_point = 3  # which point along the spline? (not min or max of u)
x_tan = central_diff(x_fine, spline_point)
y_tan = central_diff(y_fine, spline_point)
z_tan = central_diff(z_fine, spline_point)
#tangents = np.array(np.gradient([x_fine, y_fine, z_fine]))  # finds tangent direction along each point

# for single point only
# plane defined by a*x+b*y+c*z+d=0 where [a,b,c] is the normal vector.
origin = np.asarray([x_fine[spline_point], y_fine[spline_point], z_fine[spline_point]])
normal = np.array([x_tan, y_tan, z_tan])
d = np.dot(-1 * origin, normal)  # -point*normal_v

# grid points
perp_slice_size = 30  # even
precision = 0.5  # all values must be round
bounds = (perp_slice_size * precision) / 2

lowerx1 = origin[0] - bounds
upperx1 = origin[0] + bounds

lowerx2 = origin[1] - bounds

upperx2 = origin[1] + bounds

x = np.arange(lowerx1, upperx1, precision)
y = np.arange(lowerx2, upperx2, precision)

XX, YY = np.meshgrid(x, y)
ZZ = (-normal[0] * XX - normal[1] * YY - d) / normal[2];

# interpolate on CT image to generate slice
image_array = itk.GetArrayFromImage(image)
image_size = image_array.shape

interp_mesh = np.array((XX, YY, ZZ))
interp_points = np.rollaxis(interp_mesh, 0, 3)

arb_slice = scipy.interpolate.interpn((np.arange(image_size[0]), np.arange(image_size[1]), np.arange(image_size[2])),
                                      image_array, interp_points)

plt.figure()
plt.imshow(arb_slice, cmap='gray')

elapsed = time.time() - t
print('Time elapsed = ', elapsed)

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
