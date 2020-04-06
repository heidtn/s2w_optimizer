import marching_cubes
from simplefem import simplefem
import marching_cubes
from itertools import combinations 
import logging
import random
import copy
from collections import namedtuple
import numpy as np
import trimesh
from numba import njit, prange, jit
from multiprocessing import Pool

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


Parameters = namedtuple('Parameters', ['density', 'thresh', 'loads', 'constraints', 'youngs', 'poisson'])


def optimize(mesh, parameters, iterations=10):
    logger.info("Creating initial population")
    population = create_population(parameters, mesh)
    for i in range(iterations):
        logger.info("Testing population")
        fitnesses = test_population(population, parameters, mesh)
        logger.info("Peak fitness is: %s", np.max(fitnesses))
        pop_fitnesses = list(zip(population, fitnesses))
        pop_fitnesses.sort(key=lambda x: x[1], reverse=True)
        population = crossover(pop_fitnesses)
        mutate_offspring(population)

def test_population(population, parameters, mesh):
    pop_fitnesses = []
    for member in population:
        pop_fitnesses.append(single_test(member, parameters, mesh))
    return pop_fitnesses

def single_test(member, parameters, mesh):
    verts, tets, volume = member.generate_tetrahedral_hollowed_mesh(mesh, parameters.density, grid=member.grid, thresh=parameters.thresh)
    displacements, B_mats = simplefem.solve_full(tets, verts, parameters.poisson, parameters.youngs, parameters.constraints, parameters.loads)
    D = simplefem.generate_elasticity_mat(parameters.youngs, parameters.poisson)
    stresses = simplefem.get_stresses(tets, displacements, B_mats, D)
    fitness_val = fitness(stresses, volume)
    logger.info("Got new fitness: %d", fitness_val)
    return fitness_val

def create_population(parameters, mesh, pop_size=16):
    population = []
    for i in range(pop_size):
        gen = marching_cubes.MarchingCubeGenerator()
        population.append(gen)
    return population

def fitness(stresses, volume):
    # Strength to weight ratio effectively
    return np.average(np.abs(stresses)) / volume

def crossover(pop_fitnesses, parents_percentage=0.5):
    children = []
    parents = list([p[0] for p in pop_fitnesses])
    num_children = len(parents)
    parent_indices = list(combinations(range(int(round(len(parents)*parents_percentage))), 2))
    for i in range(num_children):
        parent_choices = random.choice(parent_indices)
        child = blend(parents[parent_choices[0]], parents[parent_choices[1]])
        children.append(child)
    return children

def blend(parent1, parent2, split=0.5):
    child = copy.copy(parent1)
    _fast_blend(parent1, parent2, child, split)
    return child

def _fast_blend(parent1, parent2, child, split):
    for i in range(parent1.grid.shape[0]):
        for j in range(parent1.grid.shape[1]):
            for k in range(parent1.grid.shape[2]):
                if np.random.rand() > split:
                    child.grid[i, j, k] = parent1.grid[i, j, k]
                else:
                    child.grid[i, j, k] = parent2.grid[i, j, k]

def mutate_offspring(population, prob=0.01):
    for member in population:
        mutate(member, prob)

def mutate(member, prob=0.01):
    gridsel = np.random.rand(*member.grid.shape)
    gridsel[gridsel > prob] = 0.0
    indices = np.argwhere(gridsel > 0)
    for i in indices:
        if member.grid[i[0], i[1], i[2]] > 0:
            member.grid[i[0], i[1], i[2]] = 0.0
        else:
            member.grid[i[0], i[1], i[2]] = 1.0

if __name__ == "__main__":
    constraints = [[233, [1, 1, 1]],
                   [439, [1, 1, 1]],
                   [433, [1, 1, 1]],
                   [21, [1, 1, 1]]]
    loads = [[294, 0, 0, -10.0],
             [294, 0, 0, -10.0],
             [294, 0, 0, -10.0],
             [294, 0, 0, -10.0]]
    p = Parameters(density=5, thresh=0.25, youngs=3.5, poisson=0.35, constraints=constraints, loads=loads)
    mesh = trimesh.load('featuretype.STL')
    optimize(mesh, p)