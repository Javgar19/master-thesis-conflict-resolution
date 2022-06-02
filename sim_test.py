import sys
import numpy as np
from psutil import wait_procs
import bluesky as bs
from bluesky import traffic as tr
from bluesky import settings
from bluesky.traffic.route import Route
from bluesky.navdatabase import Navdatabase
from bluesky.simulation import Simulation
from bluesky.traffic.performance.perfbase import PerfBase
import matplotlib.pyplot as plt
from bluesky.tools import geo
import random
from sacred import Experiment
import complexity_indicators as ind
import networkx as nx
from bluesky.tools.aero import nm

ex = Experiment("test-experiment")

def dist_matrix(ownship, intruder):
    I = np.eye(ownship.ntraf)

    # Horizontal conflict ------------------------------------------------------

    # qdrlst is for [i,j] qdr from i to j, from perception of ADSB and own coordinates
    qdr, dist = geo.kwikqdrdist_matrix(np.asmatrix(ownship.lat), np.asmatrix(ownship.lon),
                                np.asmatrix(intruder.lat), np.asmatrix(intruder.lon))

    # Convert back to array to allow element-wise array multiplications later on
    # Convert to meters and add large value to own/own pairs
    qdr = np.asarray(qdr)
    dist = np.asarray(dist) * nm + 1e9 * I

    return dist

def check_boundaries(traf, center, radius):
    """
    Check if any aircraft is out of the scenario bounds. It deletes it if so.
    """
    radius = radius * 1852 # From nm to meters
    id_to_delete = []
    for i in range(traf.ntraf):
        if geo.latlondist(traf.lat[i], traf.lon[i] , center[0], center[1]) > radius:
            id_to_delete.append(traf.id[i])

    if id_to_delete:
        for idx in id_to_delete:
            traf.delete(bs.traf.id.index(idx))


def init_at(radius, center, n_ac):

    for _ in range(n_ac):

        random_angle = random.random() * 360
        random_distance = radius * np.sqrt(random.random())

        random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)

        random_speed = np.random.uniform(10,20)

        angle = geo.qdrdist(random_lat, random_lon, center[0], center[1])[0]
        limit_angle = np.arccos(random_distance/(2 * radius)) * 180 / np.pi 

        acid = str(random.getrandbits(32))
        heading = random.uniform(angle - limit_angle, angle + limit_angle)

        bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon, acspd=random_speed, achdg=heading)

def init_at_wp(radius, center, n_ac, waypoints_positions):

    ac_id = 0
    for _ in range(n_ac):

        random_angle = random.random() * 360
        random_distance = radius * np.sqrt(random.random())

        random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)

        random_speed = np.random.uniform(10,20)
        #random_speed = 15

        random_wp = random.choice(waypoints_positions)
        heading = geo.qdrdist(random_lat, random_lon, random_wp[0], random_wp[1])[0]
         
        acid = str(ac_id)
        ac_id += 1

        bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon, acalt=300, acspd=random_speed, achdg=heading)
        
    return ac_id


def create_waypoints(center, radius, n_waypoints):
	"""
	The waypoints will create a polygon centered in the simulation. The function returns
	a list with each waypoint coordinates
	"""
	waypoints_positions = []

	dist_to_center =  radius/2

	for i in range(n_waypoints):
		alpha = 360/n_waypoints * i
		lat, lon = geo.qdrpos(center[0], center[1], alpha, dist_to_center)
		waypoints_positions.append([lat, lon])

	return waypoints_positions


def spawn_ac(radius, center, number_of_aircrafts):

    for _ in range(number_of_aircrafts):
        random_angle = random.random() * 360
        random_distance = radius * np.sqrt(random.random())

        random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)

        random_speed = np.random.uniform(10,20)

        angle = geo.qdrdist(random_lat, random_lon, center[0], center[1])[0]
        limit_angle = np.arccos(random_distance/(2 * radius)) * 180 / np.pi 

        acid = str(random.getrandbits(32))
        heading = random.uniform(angle - limit_angle, angle + limit_angle)

        bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon, acspd=random_speed, achdg=heading)


def spawn_ac_wp(radius, center, number_of_aircrafts, waypoints_positions, ac_id):

    for _ in range(number_of_aircrafts):
        random_angle = random.random() * 360
        random_distance = radius * np.sqrt(random.random())

        random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)

        random_speed = np.random.uniform(10,20)

        random_wp = random.choice(waypoints_positions)
        heading = geo.qdrdist(random_lat, random_lon, random_wp[0], random_wp[1])[0]
         
        acid = str(ac_id)
        ac_id += 1
        
        bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon, acalt=300, acspd=random_speed, achdg=heading)
        
    return ac_id

def front_colission(radius, center, number_of_aircrafts, waypoints_positions):

    
    random_angle = random.random() * 360
    random_distance = radius 

    random_lat1, random_lon1 = geo.qdrpos(center[0], center[1], random_angle, random_distance)
    random_lat2, random_lon2 = geo.qdrpos(center[0], center[1], random_angle + 180, random_distance)
    

    random_speed = np.random.uniform(10,20)

    random_wp = center
    heading1 = geo.qdrdist(random_lat1, random_lon1, random_wp[0], random_wp[1])[0]
    heading2 = geo.qdrdist(random_lat2, random_lon2, random_wp[0], random_wp[1])[0]
    
    ac_id = 0
    acid1 = str(ac_id)
    ac_id += 1
    acid2 = str(ac_id)
    
    bs.traf.cre(acid1, actype="M200", aclat=random_lat1, aclon=random_lon1, acalt=300, acspd=random_speed, achdg=heading1)
    bs.traf.cre(acid2, actype="M200", aclat=random_lat2, aclon=random_lon2, acalt=300, acspd=random_speed, achdg=heading2)
        
    return ac_id

def plot_at(center, radius, sources_position):

    alpha = np.linspace(0, 360, 500)


    coords = [geo.qdrpos(center[0], center[1], angle, radius) for angle in alpha]

    x = bs.traf.lat
    y = bs.traf.lon
    vx = bs.traf.gseast
    vy = bs.traf.gsnorth

    plt.figure(figsize=(8, 8))
    plt.scatter(center[0], center[1], 2)
    plt.scatter([coord[0] for coord in coords], [coord[1] for coord in coords], 1)
    plt.scatter([coord[0] for coord in sources_position], [coord[1] for coord in sources_position])
    plt.quiver(x, y, vy, vx)
    plt.title("Scenario initial configuration")
    plt.show()

def simstep():
    bs.sim.step()
    bs.net.step()
    

def change_ffmode(mode=True):
    bs.sim.ffmode = mode
    bs.sim.dtmult = 15.0

def dhij(acid1, acid2):
    """
    Horizontal distance (m) between acid1 and acid2
    """
    return geo.latlondist(bs.traf.lat[bs.traf.id.index(acid1)], bs.traf.lon[bs.traf.id.index(acid1)],
                        bs.traf.lat[bs.traf.id.index(acid2)], bs.traf.lon[bs.traf.id.index(acid2)])

def whij(acid1, acid2, H, thresh, minh):
    """
    Horizontal weight for the edge between acid1 and acid2 
    """
    dij = dhij(acid1, acid2)

    if dij <= H:
        return 1

    elif dij >= thresh:
        return 0

    else:
        return (thresh - dij)/(thresh - minh)

def tcp(acid1, acid2, Htcp, threshtcp):
    """
    Horizontal weight for the edge between acid1 and acid2 
    """
    tcpa_matrix = bs.traf.cd.tcpa
    tcpa = np.abs(tcpa_matrix[bs.traf.id.index(acid1)][bs.traf.id.index(acid2)])

    if tcpa <= Htcp:
        return 1

    elif tcpa >= threshtcp:
        return 0

    else:
        return (threshtcp - tcpa)/(threshtcp) 

def complexity_graph(H, thresh):
    graph = nx.Graph()
    graph.add_nodes_from(bs.traf.id)
    for a in graph.nodes():
        for b in graph.nodes(): # loop over every possible pair of nodes
            if a != b and (a,b) not in graph.edges():
    
                wtcp = tcp(a, b, H, thresh)
                if wtcp != 0: # Only create an edge if its weight is greater than 0
                    graph.add_edge(a, b, weight = wtcp)

    return graph

def conflict_graph():
    graph = nx.Graph()
    graph.add_nodes_from(bs.traf.id)
    for pair in bs.traf.cd.confpairs_unique:
        pair = list(pair)
        graph.add_edge(pair[0], pair[1])

    return graph

def log_variables(t, conf_graph, comp_graph, _run, num_sim, radius, n_ac, thr):
    _run.log_scalar("radius", radius)
    _run.log_scalar("n_ac", n_ac)
    _run.log_scalar("threshold", thr)
    _run.log_scalar("timestep", t)
    _run.log_scalar("num_sim", num_sim)
    _run.log_scalar("edge_density",ind.edge_density(comp_graph))
    _run.log_scalar("strength",ind.strength(comp_graph))
    #_run.log_scalar("clustering_coeff",ind.clustering_coeff(comp_graph))
    _run.log_scalar("clustering_coeff", nx.average_clustering(comp_graph, weight="weight"))
    _run.log_scalar("nn_degree",ind.nn_degree(comp_graph))
    _run.log_scalar("number_conflicts",len(bs.traf.cd.confpairs_unique))

    comp_confs = ind.comp_conf(conf_graph)
    _run.log_scalar("number_comp_conf",len(comp_confs))

    conf_lengths = [len(comp_conf) for comp_conf in comp_confs]
    if conf_lengths:
        conf_size = max(conf_lengths)
    else:
        conf_size = 0
        
    _run.log_scalar("conf_size", conf_size)



def append_variables(variables, conf_graph, comp_graph):
    """
    During a time window when conflicts are present stores the values of the variables
    in the dictionary so maximum values can be computed later
    """

    variables["ed"].append(ind.edge_density(comp_graph))
    variables["s"].append(ind.strength(comp_graph))
    variables["cc"].append(nx.average_clustering(comp_graph, weight="weight"))
    variables["nnd"].append(ind.nn_degree(comp_graph))

    for pair in bs.traf.cd.confpairs_unique:
        pair = list(pair)
        pair = {pair[0], pair[1]}
        if (not pair in variables["confs"]):
            variables["confs"].append(pair)

    comp_confs = ind.comp_conf(conf_graph)
    for conf in comp_confs:
        if not conf in variables["comp_confs"]:
            variables["comp_confs"].append(conf)

def log_conflict_variables(t, variables, _run, num_sim, radius, n_ac, thr):
    _run.log_scalar("radius", radius)
    _run.log_scalar("n_ac", n_ac)
    _run.log_scalar("threshold", thr)
    _run.log_scalar("timestep", t)
    _run.log_scalar("num_sim", num_sim)

    _run.log_scalar("edge_density", max(variables["ed"]))
    _run.log_scalar("strength",max(variables["s"]))
    _run.log_scalar("clustering_coeff",max(variables["cc"]))
    _run.log_scalar("nn_degree",max(variables["nnd"]))

    _run.log_scalar("number_conflicts", len(variables["confs"]))
    _run.log_scalar("number_comp_conf",len(variables["comp_confs"]))


    if variables["comp_confs"]:
        conf_size = max([len(comp_conf) for comp_conf in variables["comp_confs"]])
    else:
        conf_size = 0

    _run.log_scalar("conf_size", conf_size)

    if conf_size:
        print(f'Compound conflict between {conf_size} aircrafts detected')


@ex.config
def cfg():
    
    center = (47, 9)
    radius = 0.5
    n_ac = 100
    sim_time = 2*radius*1850*0.1
    n_runs = 1
    rpz = 0.13
    tcpa_thresh = 35

@ex.automain
def complexity_simulation(_run, center, radius, n_ac, sim_time, n_runs, rpz, tcpa_thresh):
    

    bs.init('sim-detached')
    ## We initialize the simulation ##

    
    #bs.traf.cd.rpz = rpz
    #bs.traf.cd.dtlookahead = 300
    
    for run in range(n_runs):
        print(f"Run {run}")
        
        bs.traf.cd.setmethod("ON")
        bs.traf.cr.setmethod("MVP")
        bs.traf.cd.rpz_def = rpz
        bs.traf.cd.dtlookahead_def = 15

        #bs.sim.simdt = 1
        #bs.sim.simt = 0
        t_max = sim_time 

        
        n_waypoints = 1

        waypoints_position = create_waypoints(center, radius, n_waypoints)
        #ac_id = init_at_wp(radius, center, n_ac, waypoints_position)
        ac_id = front_colission(radius, center, n_ac, waypoints_position)
        simstep()
        
        #plot_at(center, radius, waypoints_position)
        

        #init_at(radius, center, n_ac)

        variables = {"ed": [], "s": [], "cc": [], "nnd": [], "confs": [], "comp_confs": []}

        """ Main loop """
        change_ffmode()
        while bs.sim.simt < sim_time:

            if not bs.sim.ffmode:
                bs.sim.ffmode = True
            plot_at(center, radius, waypoints_position)
            """ Check if the acs are out of bounds and delete them if so """
            check_boundaries(bs.traf, center, radius)

            """ Spawning aircrafts """
            if bs.traf.ntraf < n_ac:
                #print(f"Before {bs.traf.ntraf} {n_ac}")
                #spawn_ac(radius, center, number_of_aircrafts = n_ac - bs.traf.ntraf)
                ac_id = spawn_ac_wp(radius, center, n_ac - bs.traf.ntraf, waypoints_position, ac_id)
                #print(f"After {bs.traf.ntraf} {n_ac}")

            conf_graph = conflict_graph()
            comp_graph = complexity_graph(4, tcpa_thresh)
           
            print(len(comp_graph.edges()))

            comp_confs = ind.comp_conf(conf_graph)
            print(len(bs.traf.cd.confpairs_unique))
            if (len(comp_confs)):
                #print(f"{len(comp_confs)} compound confs detected")
                print(comp_confs)

            if (len(bs.traf.cd.confpairs_unique) == 0) | (bs.sim.simt + bs.sim.simdt >= sim_time):

                if sum([len(variables[key]) for key in variables]) != 0:
                    log_conflict_variables(bs.sim.simt, variables, _run, run, radius, n_ac, tcpa_thresh)

                    for key in variables:
                        variables[key] = [] # reset all the values to an empty list

                
                log_variables(bs.sim.simt, conf_graph, comp_graph, _run, run, radius, n_ac, tcpa_thresh)

            else:
                append_variables(variables,  conf_graph, comp_graph)

            simstep()
            
            

        bs.sim.reset()
    
    
