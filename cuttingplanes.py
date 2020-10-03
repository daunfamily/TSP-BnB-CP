import tsplib95
import networkx as nx
from pulp import *
import time


def solve(problem, tlim):

    G = problem.get_graph()
    n = len(G.nodes)
    cuts_added = 0

    # Solution of the relaxation and set of subcycles in this solution
    cycle = []
    S = []

    # Trivial subtours with 2 vertices (i->j->i)
    for e in itertools.combinations(list(range(1, n + 1)), 1):
        S.append(e)

    # Relaxed problem (without subtour constraints)
    prob = LpProblem("TSP", LpMinimize)
    succ = LpVariable.matrix("succ", (G.nodes, G.nodes), 0, 1, LpBinary)
    # One successor and one predecessor per vertex
    for i in G.nodes:
        prob += lpSum(succ[i - 1][j - 1] for j in G.nodes) == 1
    for j in G.nodes:
        prob += lpSum(succ[i - 1][j - 1] for i in G.nodes) == 1
    prob += (
        lpSum(
            [succ[i - 1][j - 1] * G[i][j]["weight"] for j in G.nodes for i in G.nodes]
        ),
        "obj",
    )

    start = time.time()

    # Run loop while relaxed solution contains subtours
    while len(cycle) != 1:

        # Breaks if the limit of solving time is reached
        if time.time() - start > tlim:
            return ("Max time reached", None, None, tlim, None)

        # Add cuts for subtours countained in the previous solutopn of the relaxation
        for e in S:
            prob += lpSum(succ[i - 1][j - 1] for i in e for j in e) <= len(e) - 1
            cuts_added += 1
            del e

        # Solve the problem with the new constraints
        prob.solve(GUROBI(msg=0))
        print(
            "Status:",
            LpStatus[prob.status],
            "|",
            "Cuts added :",
            cuts_added,
            "|",
            "Lower bound value :",
            prob.objective.value(),
        )

        # Recover edges of the solution path in a "readable" form
        edges = []
        for i in G.nodes:
            for j in G.nodes:
                if succ[i - 1][j - 1].varValue > 0:
                    edges.append((i, j))

        # Find subtours in the solution of the relaxation
        G2 = nx.DiGraph(edges)
        cycle = [n for n in nx.simple_cycles(G2)]

        # Add the subtours on which we will add constraints
        S = []
        for i in range(len(cycle)):
            S.append([cycle[i][j] for j in range(len(cycle[i]))])

    end = time.time()

    return LpStatus[prob.status], prob.objective.value(), cycle, end - start, cuts_added


def cuttingplanes(filename, tlim):

    problem = tsplib95.load_problem("benchmarks/" + filename)

    print("-------------------")
    print("Group : Devys, Dugeon, Guyard, Hourman, Prejean")
    print("Problem :", problem.comment)
    print("-------------------")

    status, optimal_value, optimal_path, solving_time, cuts_added = solve(problem, tlim)

    print("-------------------")
    print("Problem :", problem.comment)
    print("Number of cities :", len(list(problem.get_nodes())))
    print("Method : Cutting planes")
    print("Status :", status)
    if status == "Optimal":
        print("Optimal value :", round(optimal_value))
        print("Optimal tour :", optimal_path[0])
        print("Solving time :", solving_time)
        print("Cuts added :", cuts_added)
    print("-------------------")


# --------------------------------------------------------------------------- #

if len(sys.argv) != 3:
    print("Erreur : nombre d'arguments non conforme. La commande doit être semblable à")
    print("\t python cuttingplanes.py att48.tsp 600")
else:
    if not os.path.isfile("benchmarks/" + sys.argv[1]):
        print("Le fichier de benchmarks n'existe pas")
    else:
        cuttingplanes(sys.argv[1], float(sys.argv[2]))
