import util
import numpy as np
import mc_lookups
import trimesh
import tetgen
import pyvista as pv


class VertexTracker:
    def __init__(self):
        self.bisection_index = {}
        self.vertices = []

    def get_bisect(self, vert1, vert2):
        index = (x + y for x, y in zip(vert1, vert2))
        if index in self.bisection_index.keys():
            vert_idx = self.bisection_index[index]
            return vert_idx
        else:
            self.vertices.append([(x + y)/2.0 for x, y in zip(vert1, vert2)])
            new_idx = len(self.vertices) - 1
            self.bisection_index[index] = new_idx
            return new_idx


class MarchingCubeGenerator:
    def __init__(self):
        self.grid = None
        self.vertex_tracker = VertexTracker()
        self.tris = []

    def generate_random(self, dimensions, thresh, scale_factors):
        self.grid = np.random.rand(*dimensions)
        self.dimensions = dimensions
        self.scale_factors = scale_factors
        self.grid[self.grid < thresh] = 0
        self.mask = np.ones(self.grid.shape)
        self.mask[1:-1, 1:-1, 1:-1] = 0
        self.grid += self.mask

    def clear_close_points(self, mesh, min_thresh):
        # iterate through every point, find close ones, set them to 0
        print("clearing close and external points")
        for i in range(1, self.grid.shape[0] - 1):
            for j in range(1, self.grid.shape[1] - 1):
                for k in range(1, self.grid.shape[2] - 1):
                    indices = np.array([i, j, k], dtype=np.float64)
                    indices /= np.array(self.dimensions)
                    indices -= 0.5
                    indices *= self.scale_factors
                    dist = trimesh.proximity.signed_distance(mesh, [indices])
                    if dist[0] < min_thresh:
                        self.grid[i, j, k] = 0.0
        self.mask = np.zeros(self.grid.shape)
        self.mask[1:-1, 1:-1, 1:-1] = 1.0
        self.grid *= self.mask
        print("clearing complete")


    def march(self):
        self.tris = []
        for i in range(self.grid.shape[0] - 1):
            for j in range(self.grid.shape[1] - 1):
                for k in range(self.grid.shape[2] - 1):
                    #import ipdb; ipdb.set_trace()
                    cube = self.grid[i:i+2, j:j+2, k:k+2]
                    edge_map = self.lookup_cube([i, j, k], cube)
                    if len(edge_map) >= 1:
                        for edges in edge_map:
                            vertices = []
                            for edge in edges:
                                edge0 = mc_lookups.EDGES[edge][0]
                                edge1 = mc_lookups.EDGES[edge][1]
                                vert1 = mc_lookups.VERTICES[edge0]
                                vert2 = mc_lookups.VERTICES[edge1]
                                vert1_global = [x + y for x, y in zip([i,j,k], vert1)]
                                vert2_global = [x + y for x, y in zip([i,j,k], vert2)]
                                vert_index = self.vertex_tracker.get_bisect(vert1_global, vert2_global)
                                vertices.append(vert_index)
                            self.tris.append(vertices)
        self.adjust_vertices()
    
    def adjust_vertices(self):
        vertices = np.array(self.vertex_tracker.vertices)
        vertices[:,:] /= np.array(self.dimensions, dtype=np.float64)
        vertices[:,:] -= 0.5
        vertices[:,:] *= self.scale_factors
        self.vertices = vertices
    
    def lookup_cube(self, base, vertices):
        index = sum(2**v for v in range(8) if vertices[mc_lookups.VERTICES[v]] > 0)
        edges = mc_lookups.cases[index]
        return edges



if __name__ == "__main__":
    def make_cube():
        x = np.linspace(-0.5, 0.5, 25)
        grid = pv.StructuredGrid(*np.meshgrid(x, x, x))
        return grid.extract_surface().triangulate()

    # Create to examplee PolyData meshes for boolean operations
    sphere = pv.Sphere(radius=1.65, center=(0, 0, 0))
    cube = pv.Sphere(radius=0.65, center=(0, 0, 0))

    p = pv.Plotter()
    p.add_mesh(sphere, color="yellow", opacity=0.5, show_edges=True)
    p.add_mesh(cube, color="royalblue", opacity=0.5, show_edges=True)
    p.show()

    diff = sphere.boolean_difference(cube)

    p = pv.Plotter()
    p.add_mesh(diff, opacity=0.5, show_edges=True, color=True)
    p.show()





    mesh = trimesh.load('featuretype.STL')
    mesh2 = trimesh.load('featuretype.STL')

    marcher = MarchingCubeGenerator()
    marcher.generate_random([30, 30, 30], 0.4, [6, 5, 3])
    marcher.clear_close_points(mesh, 0.15)
    marcher.march()
    
    #util.visualize_mesh(np.array(marcher.vertex_tracker.vertices), np.array(marcher.tris))
    hollow = trimesh.Trimesh(vertices=marcher.vertices, faces=marcher.tris)
    print("Hollow: ", hollow.is_watertight)

    
    manymeshes = hollow.split()
    #[trimesh.repair.fix_normals(mesh, True) for mesh in meshes]
    #print("split meshes: ", len(meshes))
    trimesh.repair.fix_normals(hollow)
    trimesh.repair.fix_normals(mesh)
    
    res = trimesh.boolean.intersection([hollow, mesh], engine='scad')
    print("res: ", res.faces)

    mesh2tet = tetgen.TetGen(mesh2.vertices, mesh2.faces)
    hollowtet = tetgen.TetGen(hollow.vertices, hollow.faces)

    print("about to subtract!")
    hollowmesh = hollowtet.mesh
    mesh2mesh = mesh2tet.mesh
    p = pv.Plotter()
    p.add_mesh(mesh2mesh, color="yellow", opacity=0.5, show_edges=True)
    p.add_mesh(hollowmesh, color="royalblue", opacity=0.5, show_edges=True)
    p.show()

    newtest = hollowmesh.boolean_union(mesh2mesh)

    #voided_part.plot()

    #hollowtet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
    #grid = hollowtet.grid
    #grid.plot(show_edges=True)

    """
    res2 = trimesh.boolean.difference([mesh2, hollow])
    print("res2: ", res2.is_watertight)

    res2.export("testout.obj")
    print("res2 type: ", type(res2.vertices))
    print("res2 verts: ", res2.vertices)

    print("Moving on!")

    tet = tetgen.TetGen(res2.vertices, res2.faces)
    dispmesh = tet.mesh
    dispmesh.plot()

    tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5, verbose=1)
    grid = tet.grid
    grid.save("testmsh.vtk")
    grid.plot(show_edges=True)    
    """

    #print(grid)
    #cells = grid.cells.reshape(-1, 5)[:, 1:]
    #cell_center = grid.points[cells].mean(1)

    """
    # extract cells below the 0 xy plane
    for i in range(10):
        mask = cell_center[:, 2] < (i+1) / 10.0 * (1.5)
        cell_ind = mask.nonzero()[0]
        subgrid = grid.extract_cells(cell_ind)

        # advanced plotting
        plotter = pv.Plotter()
        plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
        plotter.add_mesh(grid, 'r', 'wireframe')
        plotter.add_legend([[' Input Mesh ', 'r'],
                            [' Tesselated Mesh ', 'black']])
        plotter.show()
    """

    #util._add_mesh(hollow.vertices, hollow.faces)
    #util._add_mesh(mesh.vertices, mesh.faces, opacity=0.5)
    #util._add_mesh(res.vertices, res.faces)
    util._add_mesh(res2.vertices, res2.faces, opacity=0.5)
    util.mlab.show()