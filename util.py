import mayavi
from mayavi import mlab
import numpy as np
import pyvista as pv
from simplefem import simplefem
import trimesh


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

def mesh_picker(mesh_file):
    mesh = pv.read(mesh_file)
    indices = range(len(mesh.points))
    mesh["labels"] = list(indices)
    plotter = pv.Plotter()
    plotter.add_point_labels(mesh, "labels")
    plotter.add_mesh(mesh)
    plotter.show()

if __name__ == "__main__":
    mesh_picker('featuretype.STL')