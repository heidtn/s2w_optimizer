import mayavi
from mayavi import mlab
import numpy as np

def _add_mesh(vertices, faces, opacity=1.0):
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]
    
    mayavi_mesh = mlab.triangular_mesh(x, y, z, faces, opacity=opacity)

def visualize_mesh(vertices, faces, opacity=1.0):
    _add_mesh(vertices, faces, opacity)
    mlab.show()

def visualize_meshes(vertices, faces):
    for vert, face in zip(vertices, faces):
        add_mesh(vert, face)
    mlab.show()