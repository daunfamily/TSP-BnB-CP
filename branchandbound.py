import numpy as np
from copy import copy
import networkx as nx
import time
import tsplib95
import sys
import os


def compute_ub(graph):
    """ 
    Compute an upper bound based on a greddy nearest neighbour seach.
    This bound is computed at each nodes an because of the branching, the
    algorithm may fail and we return +infty in this case to force-cut the branch.
    """

    try:
        upper_bound = 0
        current_node = list(graph.nodes())[0]
        last_edges = list(graph.edges(current_node, data=True))
        while len(graph.nodes()) > 1:
            # Find nearest neighbour of current_node
            nn = sorted(
                graph.edges(current_node, data=True), key=lambda edge: edge[2]["weight"]
            )[0]
            graph.remove_node(current_node)
            upper_bound += nn[2]["weight"]
            current_node = nn[1]
        upper_bound += list(filter(lambda edge: edge[1] == current_node, last_edges))[
            0
        ][2]["weight"]
    except:
        upper_bound = float("inf")
    return upper_bound


def compute_lb(graph, ub, relax_solved, l=0):
    """
    Compute a lower bound for a tree node with the Held & Karp method (1971).
    Dual values (lambdas) can be given in parameter to begin in a "warm-state".
    """

    # M is the max number of iterations to make in the method. We set the bound at 0 first.
    M = 100
    n = len(graph.nodes())
    bound = 0

    # Lambdas of the parent node can be given in parameter to start at a warm-state to speed the resolution.
    if not l:
        l = [0 for i in range(n)]

    for m in range(M):

        # Create a subgraph without the node {1} with modified costs.
        subgraph = graph.copy()
        subgraph.remove_node(1)
        for u, v in subgraph.edges():
            subgraph[u][v]["weight"] += l[u - 1] + l[v - 1]

        # Find MST on this subgraph. If no MST found, fouce-cut the branch by setting the cost to +infty.
        try:
            tree = nx.minimum_spanning_tree(subgraph, weight="weight")
        except:
            real_cost = np.inf
            tree = None
            l = -np.inf
            break

        # Sort edges adjacent to one with modified costs.
        c = [(j, graph[1][j]["weight"] + l[0] + l[j - 1]) for j in graph.neighbors(1)]
        c = sorted(c, key=lambda v: v[1])
        e1, e2 = c[0], c[1]

        # Add two least-cost edges to the subgraph.
        tree.add_edge(1, e1[0], weight=e1[1])
        tree.add_edge(1, e2[0], weight=e2[1])

        # Compte the cost of the tree with modified edge-costs. Compute the real cost of the solution by substracting the dual values to the modified costs.
        reduced_cost = tree.size(weight="weight")
        real_cost = reduced_cost - 2 * sum(l)

        # Update the lower bound if a better bound is found.
        bound = max(bound, real_cost)

        # Update lambda according to Volgenant and Jonker's method.
        # (https://www.sciencedirect.com/science/article/pii/S0377221796002147)
        t1 = real_cost / (2 * n)
        t = (
            t1 * (m - 1) * (2 * M - 5) / (2 * M - 2)
            - t1 * (m - 2)
            + 0.5 * t1 * (m - 1) * (m - 2) / ((M - 1) * (M - 2))
        )
        new_l = [l[i] + t * (tree.degree(i + 1) - 2) for i in range(len(l))]

        # If the lambdas don't change a lot, we stop the relaxation.
        if np.linalg.norm(np.array(l) - np.array(new_l), np.inf) < 10 ** (-4):
            break
        l = new_l
        relax_solved += 1

    return (real_cost, tree, l, relax_solved)


def cycle_detection(tree):
    """
    Return wether the 1-tree is a cycle or not.
    """

    try:
        return len(nx.find_cycle(tree)) == len(tree.nodes())
    except:
        return False


def find_branching(path):
    """
    Find a vertex with more than 2 incident edges and return its inscident edges.
    """

    for node in path.nodes():
        if path.degree(node) > 2:
            for edge in nx.edges(path, node):
                yield edge
            return None


def solve(problem, tlim):
    """
    Solve TSP with 1-tree Lagrangian relaxation (Held & Karp method to compute a lower bound).
    Each node of the tree has for parameter :
        - Its dept in the tree
        - A graph derived form the initial graph (with some edges removed)
        - The path and the cost of the computed lower bound 
        - Dual values (lambdas) of the relaxation which are used to start the lower bound computation in a "warm state"
    """

    # Tsplib95 add a loop over each node by default so we ave to remove all these loops.
    graph = problem.get_graph()
    for node in problem.get_nodes():
        graph.remove_edge(node, node)

    # Relabel nodes which should starts at 1.
    mapping = dict(zip(graph.nodes(), range(1, len(graph.nodes()) + 1)))
    graph = nx.relabel_nodes(graph, mapping)

    start = time.time()

    ### Initialisation of the tree ###

    # Store the nodes to be explored in a pile. The first upper bound is computed with a greedy nearest-neighbour algorithm.
    node_to_explore = []
    ub = compute_ub(graph.copy())
    print("First ub with nearest neighbour algorithm :", ub)
    best_path = None
    nodes_explored = 0
    relax_solved = 0

    # Create the root with the initial graph with all lambda set to 0.
    # Compute a lower bound and test if the solution of the relaxation is a feasible solution better than the actual solution.
    root = {"graph": graph, "level": 0}
    root["lb"], root["path"], root["l"], relax_solved = compute_lb(
        root["graph"].copy(), ub, relax_solved
    )
    if cycle_detection(root["path"]) and root["lb"] < ub:
        ub = root["lb"]
        best_path = root["path"]

    print(
        "Nodes explored : 0",
        "| Resolved relaxations :",
        relax_solved,
        "| Upper bound :",
        format(ub, ".2f"),
        "| Lower bound :",
        format(root["lb"], ".2f"),
        "| Gap :",
        format(100 * abs(ub - root["lb"]) / ub, ".2f"),
        "%",
    )

    # Initialize the exploration pile.
    node_to_explore.append(root)

    ### Exploration of the tree ###

    while node_to_explore:

        # Breaks if the limit of solving time is reached
        if time.time() - start > tlim:
            return ("Max time reached", None, None, tlim, None, None)

        # Take the node on the top of the pile
        node = node_to_explore.pop()
        nodes_explored += 1
        if nodes_explored % 10 == 0:
            print(
                "Nodes explored :",
                nodes_explored,
                "| Resolved relaxations :",
                relax_solved,
                "| Upper bound :",
                format(ub, ".2f"),
                "| Current lower bound :",
                format(node["lb"], ".2f"),
                "| Gap :",
                format(100 * abs(ub - node["lb"]) / ub, ".2f"),
                "%",
            )

        # Find the vertex to branch on and the edges to remove
        edges_to_remove = list(find_branching(node["path"]))

        # Create one child per edges incident to this vertex
        child_nodes = [{} for edge in edges_to_remove]

        for edge in range(len(edges_to_remove)):

            # Create a child but remove the edge corresponding to this child from the graph.
            child_nodes[edge] = {
                "graph": node["graph"].copy(),
                "level": node["level"] + 1,
            }
            child_nodes[edge]["graph"].remove_edge(
                edges_to_remove[edge][0], edges_to_remove[edge][1]
            )

            # Test if the graph without the edge is still connected.
            if nx.is_connected(child_nodes[edge]["graph"]):

                # Compute an upper bound with neirest-neighbour heurisitic and try to update the global upper-bound.
                local_ub = compute_ub(child_nodes[edge]["graph"].copy())
                ub = min(ub, local_ub)

                # Compute a lower bound and the path associated to this bound according to the Held & Karp method.
                # Gives the lambdas of the parent node to start the relaxation in a "warm state".
                (
                    child_nodes[edge]["lb"],
                    child_nodes[edge]["path"],
                    child_nodes[edge]["l"],
                    relax_solved,
                ) = compute_lb(
                    child_nodes[edge]["graph"].copy(), local_ub, relax_solved, node["l"]
                )

                # Test if it's a new upper bound.
                if (
                    cycle_detection(child_nodes[edge]["path"])
                    and child_nodes[edge]["lb"] < ub
                ):
                    ub = child_nodes[edge]["lb"]
                    best_path = child_nodes[edge]["path"]

            else:
                # If graph is not connected, problem is infeasible for the current node so we force-cut the branch.
                child_nodes[edge]["lb"] = np.inf

        # Keep only childs with lb < ub and add with best-first strategy.
        child_nodes = [
            child_node for child_node in child_nodes if child_node["lb"] < ub
        ]
        node_to_explore.extend(sorted(child_nodes, key=lambda node: node["lb"]))

    stop = time.time()

    if best_path:
        status = "Optimal"
        cycle = nx.find_cycle(best_path)
    else:
        status = "Unsolved"

    return status, ub, cycle, stop - start, nodes_explored, relax_solved


def branchandbound(filename, tlim):

    problem = tsplib95.load_problem("benchmarks/" + filename)

    print("-------------------")
    print("Group : Devys, Dugeon, Guyard, Hourman, Prejean")
    print("Problem :", problem.comment)
    print("-------------------")

    (
        status,
        optimal_value,
        optimal_path,
        solving_time,
        nodes_explored,
        relax_solved,
    ) = solve(problem, tlim)

    print("-------------------")
    print("Problem :", problem.comment)
    print("Number of cities :", len(list(problem.get_nodes())))
    print("Method : Branch and Bound")
    print("Status :", status)
    if status == "Optimal":
        print("Optimal value :", round(optimal_value))
        print("Optimal tour :", optimal_path)
        print("Solving time :", solving_time)
        print("Nodes explored :", nodes_explored)
        print("Resolved relaxations :", relax_solved)
    print("-------------------")


# ------------------------------------------------------------- #

if len(sys.argv) != 3:
    print("Erreur : nombre d'arguments non conforme. La commande doit être semblable à")
    print("\t python branchandbound.py att48.tsp 600")
else:
    if not os.path.isfile("benchmarks/" + sys.argv[1]):
        print("Le fichier de benchmarks n'existe pas")
    else:
        branchandbound(sys.argv[1], float(sys.argv[2]))
