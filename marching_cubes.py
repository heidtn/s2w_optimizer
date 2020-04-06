import util
import numpy as np
import mc_lookups
import trimesh
import tetgen
import pyvista as pv
from simplefem import simplefem
import pygmsh


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

    def __copy__(self):
        new = MarchingCubeGenerator()
        new.set_grid(self.grid, self.thresh, self.scale_factors, self.density, self.center)
        return new

    def generate_random(self, thresh, extents, density, center):
        count = np.round(np.array(extents) * density).astype(np.int)
        self.grid = np.random.rand(*count)
        self.set_grid(self.grid, thresh, extents, density, center)

    def set_grid(self, grid, thresh, extents, density, center):
        self.thresh = thresh
        self.grid = np.array(grid)
        count = np.round(np.array(extents) * density).astype(np.int)
        self.dimensions = count
        self.scale_factors = extents
        self.grid[self.grid < thresh] = 0
        self.mask = np.ones(self.grid.shape)
        self.mask[1:-1, 1:-1, 1:-1] = 0
        self.grid += self.mask
        self.center = center
        self.density = density

    def clear_close_points(self, mesh, min_thresh):
        # iterate through every point, find close ones, set them to 0
        print("clearing close and external points")
        for i in range(1, self.grid.shape[0] - 1):
            for j in range(1, self.grid.shape[1] - 1):
                for k in range(1, self.grid.shape[2] - 1):
                    indices = np.array([i, j, k], dtype=np.float64)
                    indices /= self.dimensions
                    indices -= 0.5
                    indices *= self.scale_factors
                    indices += self.center
                    dist = trimesh.proximity.signed_distance(mesh, [indices])
                    if dist[0] < min_thresh:
                        self.grid[i, j, k] = 0.0
        self.mask = np.zeros(self.grid.shape)
        self.mask[1:-1, 1:-1, 1:-1] = 1.0
        self.grid *= self.mask
        print("clearing complete: ")


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

    def generate_tetrahedral_hollowed_mesh(self, mesh, density, grid=None, thresh=0.25):
        # density in points per square unit, 
        # Does the whole generation process. If grid is None, generates a random one
        extents = mesh.bounding_box.extents
        center = mesh.bounds.mean(axis=0)
        if grid is not None:
            self.set_grid(grid, thresh, extents, density, center)
        else:
            self.generate_random(thresh, extents, density, center)
        self.clear_close_points(mesh, thresh)
        self.march()

        hollow = trimesh.Trimesh(vertices=self.vertices, faces=self.tris)

        verts, tets = generate_dual_tetrahedrons(mesh, hollow)

        new_volume = mesh.volume - hollow.volume
        return verts, tets, new_volume
    
    def adjust_vertices(self):
        vertices = np.array(self.vertex_tracker.vertices)
        vertices[:,:] /= np.array(self.dimensions, dtype=np.float64)
        vertices[:,:] -= 0.5
        vertices[:,:] *= self.scale_factors
        vertices[:,:] += self.center
        self.vertices = vertices
    
    def lookup_cube(self, base, vertices):
        index = sum(2**v for v in range(8) if vertices[mc_lookups.VERTICES[v]] > 0)
        edges = mc_lookups.cases[index]
        return edges


def remove_tets_internal(internal_mesh, tets, vertices):
    distances = trimesh.proximity.signed_distance(internal_mesh, vertices)
    indices = np.argwhere(distances > 0)
    index_set = set(indices.T.tolist()[0])
    new_tets = []
    for tet in tets:
        valid_node = True
        for node in tet:
            if node in index_set:
                valid_node = False
        if valid_node:
            new_tets.append(tet)
    print("removed: ", len(tets) - len(new_tets), " elements")
    return np.array(new_tets)

def remove_tets_convex(internal_mesh, tets, vertices):
    # Need to remove tets that are inside the internal mesh, but
    # not tests that are connected outside.  So any tet whose faces are all faces in the 
    # Internal mesh should be removed.
    new_tets = []
    mesh_tris = set([frozenset(tri) for tri in internal_mesh.faces])
    for tet in tets:
        tet_wrap = np.append(tet, tet)
        tris = set([frozenset(tet_wrap[i:i+3]) for i in range(4)])
        if not tris.issubset(mesh_tris):
            new_tets.append(tet)
    return new_tets

def generate_dual_tetrahedrons(mesh, internal_mesh):
    # Meshes are of trimesh type
    diff = trimesh.boolean.difference((mesh, internal_mesh))
    new_verts = diff.vertices
    new_faces = diff.faces

    tet = tetgen.TetGen(new_verts, new_faces)
    tet.tetrahedralize(order=1)
    grid = tet.grid
    tets = simplefem.extract_tets(grid.cells)
    verts = grid.points

    #new_tets = remove_tets_internal(internal_mesh, tets, verts)
    new_tets = remove_tets_convex(internal_mesh, tets, verts)
    verts = np.array(verts)
    new_tets = np.array(new_tets)
    #simplefem.display_tets(verts, new_tets)
    return verts, new_tets

def test_full_gen():
    bigcube = trimesh.creation.box((1, 1, 1))
    smallcube = trimesh.creation.box((0.5, 0.5, 0.5))
    generate_dual_tetrahedrons(bigcube, smallcube)

def test_marcher():
    mesh = trimesh.load('featuretype.STL')
    extents = mesh.bounding_box.extents
    center = mesh.bounds.mean(axis=0)
    print("extents: ", extents)
    marcher = MarchingCubeGenerator()
    marcher.generate_random(0.25, extents, 10.0, center)
    marcher.clear_close_points(mesh, 0.15)
    marcher.march()
    hollow = trimesh.Trimesh(vertices=marcher.vertices, faces=marcher.tris)
    hollow.show()
    diff = trimesh.boolean.difference((mesh, hollow))
    print("diff volume: ", diff.volume)
    print("mesh volume: ", mesh.volume)
    diff.show()
    util.visualize_mesh(diff.vertices, diff.faces)

if __name__ == "__main__":
    test_marcher()