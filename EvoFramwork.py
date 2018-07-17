import random
import numpy
import time
import evaluation as ev
import networkx as nx

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

# Constants
init_pop_size = 300
number_of_nodes = ev.GRAPH_SIZE
# squared because that many elements in adj matrix
ind_size = number_of_nodes * number_of_nodes
# CXPB  is the probability with which two individuals
# are crossed
CXPB = 0.5
# MUTPB is the probability for mutating an individual
MUTPB = 0.4
# NGEN  is the number of generations for which the
# evolution runs
NGEN =100000

# Deap ---------------------------------

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

# Attribute generator
# define 'attr_bool' to be an attribute ('gene')
# which corresponds to integers sampled uniformly
# from the range [0,1] (i.e. 0 or 1 with equal
# probability)
toolbox.register("attr_bool", random.randint, 0, 1)

# 'Individual' is an individual
# consisting of 100 'attr_bool' elements ('genes')
toolbox.register("individual", tools.initRepeat,
                 creator.Individual, toolbox.attr_bool, n=ind_size)

# The population is a list of individuals
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


# ---------------------------------------------------------


def cxTwoPointCopy(ind1, ind2):
    """Execute a two points crossover with copy on the input individuals. The
    copy is required because the slicing in numpy returns a view of the data,
    which leads to a self overwritting in the swap operation. It prevents
    ::

        import numpy
        a = numpy.array((1,2,3,4))
        b = numpy.array((5.6.7.8))
        a[1:3], b[1:3] = b[1:3], a[1:3]
        print(a)
        >[1 6 7 4]
        print(b)
        >[5 6 7 8]
    """
    size = len(ind1)
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else:  # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()

    return ind1, ind2


# ----------
# Operator registration
# ----------
# register the goal / fitness function
toolbox.register("evaluate", ev.evaluate)

# register the crossover operator
toolbox.register("mate", tools.cxTwoPoint)

# register a mutation operator with a probability to
# flip each attribute/gene of 0.05
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)

# operator for selecting individuals for breeding the next
# generation: each individual of the current generation
# is replaced by the 'fittest' (best) of three individuals
# drawn randomly from the current generation.
toolbox.register("select", tools.selTournament, tournsize=3)


# Make sure that the trace is 0.
# Manually change everything in diagonal to 0.
def fix_pop(pop):
    for p in pop:
        # change the correct value in the list.
        # Start at zero and change every n+1'th one,
        # where n is the number of nodes
        for i in range(0, ind_size, ev.GRAPH_SIZE + 1):
            p[i] = 0


# Takes in a list and outputs gephi file


def output_graph(ind):
    # a is the adjacency matrix
    a = numpy.reshape(ind, (ev.GRAPH_SIZE, ev.GRAPH_SIZE))
    D = nx.DiGraph(a)
    nx.write_graphml(
        D, "{} with {} iterations.graphml".format(ev.GRAPH_SIZE, NGEN))
    print("Graph Saved.")


# ---------------------------------------------------------------


def main():
    # random.seed(64)
    print("Running Evolution...")

    pop = toolbox.population(n=init_pop_size)

    # Make pop have trace of 0
    fix_pop(pop)

    # Numpy equality function (operators.eq) between two arrays returns the
    # equality element wise, which raises an exception in the if similar()
    # check of the hall of fame. Using a different equality function like
    # numpy.array_equal or numpy.allclose solve this issue.

    hof = tools.HallOfFame(10, similar=numpy.array_equal)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    # evolve(pop, stats)

    # CXPB  is the probability with which two individuals
    #       are crossed
    #
    # MUTPB is the probability for mutating an individual
    #
    # NGEN  is the number of generations for which the
    #       evolution runs
    algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB,
                        ngen=NGEN, stats=stats, halloffame=hof, verbose=True)

    return pop, stats, hof


if __name__ == "__main__":
    show_generation = False
    show_pop_size = False
    show_stats = False

    print("Evolution start with initial population size {} and individual size {}".format(
        init_pop_size, ev.GRAPH_SIZE))

    start_time = time.time()
    data = main()
    print("--- %s seconds ---" % (time.time() - start_time))

    # Make vars for results
    current_pop = data[0]
    stats = data[1]
    hall = data[2]

    # Results
    print("\nResults:")
    best_ind = tools.selBest(data[0], 1)[0]
    print("Best individual is \n %s, %s" %
          (ev.adj_matrix(best_ind), best_ind.fitness.values))
    ev.write_best(best_ind)

    # Output graph to Gephi
    print("")
    output_graph(list(best_ind))


# ----------- Old Functions ---------------#

# -------------------- Custom Evo --------------------------------

def evolve(pop, stats):
    # CXPB  is the probability with which two individuals
    #       are crossed
    #
    # MUTPB is the probability for mutating an individual
    #
    # NGEN  is the number of generations for which the
    #       evolution runs
    CXPB, MUTPB, NGEN = 0.5, 0.2, 40

    print("Start of evolution")

    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    if show_pop_size:
        print("  Evaluated %i individuals" % len(pop))

    # Begin the evolution
    for g in range(NGEN):
        if show_generation:
            print("-- Generation %i --" % g)

        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))

        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):

            # cross two individuals with probability CXPB
            if random.random() < CXPB:
                toolbox.mate(child1, child2)

                # fitness values of the children
                # must be recalculated later
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:

            # mutate an individual with probability MUTPB
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        if show_pop_size:
            print("  Evaluated %i individuals" % len(invalid_ind))

        # The population is entirely replaced by the offspring
        pop[:] = offspring

        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]

        if show_stats:
            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x * x for x in fits)
            std = abs(sum2 / length - mean ** 2) ** 0.5

            print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            print("  Std %s" % std)

        record = stats.compile(pop)
        print(record)

    print("-- End of (successful) evolution --")

    return pop
