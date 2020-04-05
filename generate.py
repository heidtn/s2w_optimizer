import marching_cubes.py
from simplefem import simplefem


class GeneticMember:
    def __init__(self):
        pass


def optimize(mesh, parameters, iterations=50):
    population = create_population(parameters)
    for i in range(iterations):
        results = test_population(population, parameters, mesh)
        # log the results
        population = mate_population(population, results, parameters, mesh)


def test_population(population, parameters, mesh):
    pass

def mate_population(population, results, parameters, mesh):
    pass

def create_population(parameters, pop_size=16):
    pass

def initialize(mesh, parameters):
    pass

def fitness(strains):
    pass

def crossover(gene1, gene2):
    pass

def mutate(gene):
    pass

def create_next_generation(population):
    pass

if __name__ == "__main__":
    # take in a mesh and load config
    # maybe a density config?
    # iterate on that bish