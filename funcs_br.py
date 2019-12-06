import glob
import itk
import os
import sys
import time

import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import interpolate
from scipy import ndimage as nd
from skimage import measure, feature
from skimage.morphology import skeletonize, medial_axis


def save_image(array, affine, header, namepath):
    # This is only for images loaded with nibabel or numpy arrays
    nib.save(nib.Nifti1Image(array, affine, header), namepath + '.nii.gz')
    return

def plot_3d(*args, threshold=0):
    ### POSITIONAL ARGS respectively:
    # required: image.
    # optional: graph verticies
    
    # TODO - figure the proper orientation out
    # Position the scan upright,
    # so the head of the patient would be at the top facing the camera
    p = args[0].T #.transpose(2, 1, 0)
    # p = p[:, :, ::-1]

    verts, faces = measure.marching_cubes_classic(p, threshold)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces], alpha=0.1)
    face_color = [0.5, 0.5, 1]
    mesh.set_facecolor(face_color)
    ax.add_collection3d(mesh)
    if len(args) == 2:
        graph_vertices = args[1]
        ax.scatter(graph_vertices[..., 2], graph_vertices[..., 1], graph_vertices[..., 0], s=30, c='red', marker='o')

    ax.set_xlim(0, p.shape[0])
    ax.set_ylim(0, p.shape[1])
    ax.set_zlim(0, p.shape[2])

    plt.show()
    
    
def isosurface(binary_array):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    pos = np.where(binary_array==1)
    
    ax.scatter(pos[0], pos[1], pos[2], c='black')
    ax.set_xlim(0, binary_array.shape[0])
    ax.set_ylim(0, binary_array.shape[1])
    ax.set_zlim(0, binary_array.shape[2])
    plt.show()
