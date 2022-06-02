import networkx as nx
import rand_conflict

def comp_conf(conf_pairs):
	"""
	Return a list with all compound conficts
	"""
	comp_conf = [{conf_pairs[0][0], conf_pairs[0][1]}]
	for pair in conf_pairs:
		n = 0

		for i, conflict in enumerate(comp_conf):
			
			if (pair[0] not in conflict) & (pair[1] not in conflict):
				n += 1
			else:
				conflict.add(pair[0])
				conflict.add(pair[1])
				break

			if n == len(comp_conf):
				comp_conf.append({pair[0], pair[1]})

				
	return comp_conf

#print(comp_conf([[1,2], [3,4], [5, 3], [2, 6], [3, 2]]))

def comp_conf_nx(acids, conf_pairs):

	g = nx.Graph()
	g.add_nodes_from(acids)

	for pair in conf_pairs:
		g.add_edge(pair[0], pair[1])

	d = list(nx.connected_components(g)) 
	return d


r = rand_conflict.generate(4, 2000, 3000, 15, 8)
rand_conflict.plot_at(r)

