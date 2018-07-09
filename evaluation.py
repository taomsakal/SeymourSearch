import sys
import numpy as np

# Constants

GRAPH_SIZE = 18


# Evaluate the individual
# That is, get it's fitness
# The higher the average anti-satisfaction, the higher the fitness.
def evaluate(ind):
    # Modify the actual individual to have trace of zero
    fix_individual(ind)

    individual = np.copy(ind)  # Copy so no errors

    num_edges = sum(individual)

    # Get the adjacency matrix
    A = adj_matrix(individual)

    # if has self loops, then zero fitness
    # and immediatly return
    if self_loops(A):
        fitness = 0
        return fitness,  # Comma is important

    A_sq = np.linalg.matrix_power(A, 2)

    # Loose fitness for having diagons
    # Don't set to zero or else almost all will have diagons
    # and pop will never evolve.
    fit_minus = diagons(A_sq)

    # Get the nhood matrix
    nhood = nhood_matrix(A, A_sq)

    # Get the vector of antisatisfactions
    sv = sum_vector(nhood)

    # See if individual is perfect. (No diagons)
    # If so, abort program, return individual, and write individual to file.
    if fit_minus == 0 and is_perfect(sv):
        success(A)

    # Make the fitness the average anti_satisfaction
    # fitness = avg_anti_satisfaction(sv) - fit_minus
    # fitness = avg_anti_satisfaction(sv) - fit_minus + (tt_number(A))*0.01
    # fitness = (avg_anti_satisfaction(sv)*GRAPH_SIZE - fit_minus + tt_number(A)*GRAPH_SIZE
    #        - sum(range(0, num_edges)) - sink_number(A)*GRAPH_SIZE - num_sym_points(sv)*GRAPH_SIZE*GRAPH_SIZE)

    fitness = symr_diamond_number(
        A_sq) - fit_minus + avg_anti_satisfaction(sv) * GRAPH_SIZE - sum(range(0, num_edges)) * .001

    return fitness,  # <-- The comma is important!


def fix_individual(ind):
    for i in range(0, GRAPH_SIZE * GRAPH_SIZE, GRAPH_SIZE + 1):
        ind[i] = 0


def triangle_number(A):
    A_cub = np.linalg.matrix_power(A, 3)
    return np.trace(A_cub) / 3


def tt_number(A):
    A_cub = np.linalg.matrix_power(A, 3)
    A_trans = np.transpose(A)
    ntrans = np.trace(np.linalg.matrix_power(
        A + A_trans, 3)) / 6 - (np.trace(A_cub) / 3)
    return ntrans


# Get the number of seymour diamonds.
# These diamonds are two paths of length two going from node i to j.
# So to find the number of them we take all entries that are bigger than two
# (these have at least two paths, so at least one diamond) and out of that
# number choose two to get the number of diamonds from that vertex.
# Then sum up the vertices.


def symr_diamond_number(A_sq):
    count = 0
    for r in A_sq:
        for i in r:
            if i > 1:
                count += choose(i, 2)

    return count


def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):  # changed from xrange
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


def sink_number(A):
    out_degree = sum_vector(A)
    sink_num = 0
    for v in out_degree:
        if v == 0:
            sink_num += 1

    return sink_num


# Get the number of seymour points


def num_sym_points(sv):
    count = 0
    for v in sv:
        if v <= 0:
            count += 1

    return count


def min_ant_sat(sv):
    min = 1000000
    for v in sv:
        if v < min:
            min = v
    return min


def success(individual):
    write(individual)  # Make text file with the perfect one
    sys.exit("We have found a perfect individual. It is\n" + str(individual))


def write(individual):
    f = open('Perfection', 'w')
    f.write(str(individual))  # python will convert \n to os.linesep
    f.close()  # you can omit in most cases as the destructor will call it


def write_best(individual):
    f = open("Best for" + str(GRAPH_SIZE), 'w')
    f.write(str(individual))  # python will convert \n to os.linesep
    f.close()  # you can omit in most cases as the destructor will call it


# Make the adjacency matrix from the


def adj_matrix(individual):
    # Reshape into square numpy array
    return np.reshape(individual, (GRAPH_SIZE, GRAPH_SIZE))


# We make the neighborhood matrix, with values of 1 at first neighbors and
# values of -1 for the second neighborhoods.
# A is the adjacency matrix


def nhood_matrix(A, A_squared):
    Nhood = np.copy(A)

    # This is the overlapping step
    for i in range(0, GRAPH_SIZE):
        for j in range(0, GRAPH_SIZE):
            # If path of length 2 and no path of length 1
            # place a -1 to show is in second nhood
            if A_squared[i][j] != 0 and A[i][j] != 1:
                Nhood[i][j] = -1
        # This takes care of all the cases. Any others either have no
        # path or already has a path of length 1, which Nhood already has
        # since it is set to A

    return Nhood


# Vector of the antisatisfaction of each node, calculated via the nhood matrix


def sum_vector(nhood_matrix):
    sum_vector = []
    for i in range(0, GRAPH_SIZE):
        sum_vector.append(np.sum(nhood_matrix[i]))

    return sum_vector


# Make our program abort if it actually does find a non-seymour graph.


def is_perfect(sum_vector):
    perfect = True
    # Go through each antisatisfaction. If any is
    # negative or 0, then that is a seymour point and
    # we set perfect to false.
    for i in sum_vector:
        if i <= 0:
            perfect = False

    if perfect:
        return True
    return False


# tests if there are self loops


def self_loops(A):
    # Since all elements are positive, just test the trace
    if np.trace(A) == 0:
        return False

    return True


def has_sink(A):
    # If any row is all 0, then has sink
    sv = sum_vector(A)


# Tests for diagons
# Returns the number of them


def diagons(A_sq):
    # No paths of length 2 to self if no diagons
    return np.trace(A_sq) * GRAPH_SIZE * GRAPH_SIZE  # Huge selection against


# Return our average anti_satisfaction


def avg_anti_satisfaction(sum_vector):
    return float(sum(sum_vector)) / float(len(sum_vector))


def test(individual):
    fitness = 0

    # Get the adjacency matrix
    A = adj_matrix(individual)
    print("THis is matrix A")
    print(A)

    # if has self loops, then zero fitness
    # and immediatly return
    if self_loops(A):
        print("self loops")
        fitness = 0
        return fitness,  # Comma is important

    print("No self loops")

    A_sq = np.dot(A, A)

    # if diagons(A_sq):
    #    print("diagons")
    #    fitness = 0
    #    return fitness

    # Get the nhood matrix
    Nhood = nhood_matrix(A, A_sq)
    print(Nhood)

    # Get the vector of antisatisfactions
    sv = sum_vector(Nhood)
    print(sv)

    # See if individual is perfect.
    # If so, abort program, return individual, and write individual to file.
    if is_perfect(sv):
        print("perfect!")
        success(individual)

    # Make the fitness the average anti_satisfaction
    fitness = avg_anti_satisfaction(sv)

    print(fitness)

    print("A returns")
    print(A)
    return fitness,  # <-- The comma is important!
