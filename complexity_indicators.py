import numpy as np
import networkx as nx


def edge_density(graph):

    # Sum of the weights 
    sum_w = sum([d["weight"] for _, _, d in G.edges(data=True)])
    V = G.number_of_nodes()
    return sum_w / (V*(V-1)/2)

def strength(graph, mean=True):
	"""
	Return the strength of every node or the mean value if 'mean' parameter is True
	"""

	# Weighted node degrees
	node_degrees = graph.degree(weight = "weight")

	if mean is False:
		return node_degrees

	# Mean of the degrees 
	degrees_mean = np.mean([degree[1] for degree in node_degrees])
	return degrees_mean

def node_triplets(G, node):
	"""
	Auxiliary funtion that computes and returns the triplets that 'node' belongs to
	"""
	triplets = []
	for neighbor in G.adj[node]:
		for adj in G.adj[neighbor]:
			if adj in G.adj[node]:
				triplet = {node, neighbor, adj}
				if triplet not in triplets:
					triplets.append(triplet)
	return [list(triplet) for triplet in triplets]

def clustering_coeff(G, mean=True):
    """
    Return the clustering coefficient for every node or the mean value if 'mean' parameter is True
    """
    clustering_coeffs = []
    for node in G.nodes:
        triplets = node_triplets(G, node) # Compute all triplets that node is part of
        
        if triplets:
            
            wijk = 0
            for triplet in triplets:
                wijk += G[triplet.pop(triplet.index(node))][triplet[0]]["weight"] + G[triplet[0]][triplet[1]]["weight"]

            clustering_coeffs.append(wijk/(2 * G.degree(weight = "weight", nbunch=node) * (G.degree(nbunch=node) - 1)))

        else:
            clustering_coeffs.append(0)
            
    if mean is False:
        return clustering_coeffs
    
    return np.mean(clustering_coeffs)


def nn_degree(graph):
	nodes_nn_degree = list(nx.average_neighbor_degree(G, weight = "weight").values())
	return np.mean(nodes_nn_degree)
	


