import util
import numpy as np
import mc_lookups
import trimesh


class VertexTracker:
    def __init__(self, scale, offset):
        self.bisection_index = {}
        self.vertices = []
        self.scale = scale
        self.offset = offset

    def get_bisect(self, vert1, vert2):
        index = (x + y for x, y in zip(vert1, vert2))
        if index in self.bisection_index.keys():
            vert_idx = self.bisection_index[index]
            return vert_idx
        else:
            vert = (np.array(vert1) + np.array(vert2)) / 2.0
            vert *= self.scale
            vert += self.offset
            self.vertices.append(vert.tolist())
            new_idx = len(self.vertices) - 1
            self.bisection_index[index] = new_idx
            return new_idx


class MarchingCubeGenerator:
    def __init__(self, scale_factors, offset):
        self.grid = None
        self.scale_factors = scale_factors
        self.vertex_tracker = VertexTracker(scale_factors, offset)
        self.tris = []

    def generate_random(self, dimensions, thresh):
        self.grid = np.random.rand(*dimensions)
        self.dimensions = dimensions
        self.grid[self.grid < thresh] = 0
        self.mask = np.ones(self.grid.shape)
        self.mask[1:-1, 1:-1, 1:-1] = 0
        self.grid += self.mask

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

    def elimate_close_points(self, mesh):




if __name__ == "__main__":
    marcher = MarchingCubeGenerator()
    marcher.generate_random([30, 30, 30], 0.4, [6, 5, 3])
    marcher.march()
    
    #util.visualize_mesh(np.array(marcher.vertex_tracker.vertices), np.array(marcher.tris))
    mesh = trimesh.load('featuretype.STL')
    mesh2 = trimesh.load('featuretype.STL')
    hollow = trimesh.Trimesh(vertices=marcher.vertices, faces=marcher.tris)
    print("Hollow: ", hollow.is_watertight)

    
    #meshes = hollow.split()
    #[trimesh.repair.fix_normals(mesh, True) for mesh in meshes]
    #print("split meshes: ", len(meshes))
    trimesh.repair.fix_normals(hollow)
    trimesh.repair.fix_normals(mesh)
    res = trimesh.boolean.intersection([hollow, mesh], engine='scad')
    print("res: ", res.faces)

    res2 = trimesh.boolean.difference([mesh2, res])
    print("res2: ", res2.faces)

    #util._add_mesh(mesh.vertices, mesh.faces)
    #util._add_mesh(res.vertices, res.faces)
    util._add_mesh(res2.vertices, res2.faces)
    util.mlab.show()