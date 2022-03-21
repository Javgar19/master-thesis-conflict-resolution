import sys
import numpy as np
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

ex = Experiment("test-experiment")


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


def create_sources(center, radius, n_sources):
    """
    The sources will create a polygon centered in the simulation. The function returns
    a list with each source coordinates
    """
    sources_positions = []

    dist_to_center = radius

    for i in range(n_sources):
        alpha = 360/n_sources * i
        lat, lon = geo.qdrpos(center[0], center[1], alpha, dist_to_center)
        sources_positions.append([lat, lon])

    return sources_positions



def init_at(radius, center, n_ac):

    for _ in range(n_ac):

        random_angle = random.random() * 360
        random_distance = radius * np.sqrt(random.random())

        random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)


        angle = geo.qdrdist(random_lat, random_lon, center[0], center[1])[0]
        limit_angle = np.arccos(random_distance/(2 * radius)) * 180 / np.pi 

        acid = str(random.getrandbits(32))
        heading = random.uniform(angle - limit_angle, angle + limit_angle)

        bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon, acspd=40, achdg=heading)


def spawn_ac(sources_position, radius, center, number_of_aircrafts):

    for _ in range(number_of_aircrafts):
        source = random.choice(sources_position)
        angle = geo.qdrdist(source[0], source[1], center[0], center[1])[0]
        limit_angle = 30

        acid = str(random.getrandbits(32))
        heading = random.uniform(angle - limit_angle, angle + limit_angle)

        bs.traf.cre(acid, actype="M200", aclat=source[0], aclon=source[1], acspd=40, achdg=heading)


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
        bs.sim.dtmult = 5.0

def at_to_graph():
    graph = nx.Graph()
    graph.add_nodes_from(bs.traf.id)
    for pair in bs.traf.cd.confpairs_all:
        graph.add_edge(pair[0], pair[1])

    return graph

def log_variables(t, graph, _run):
    _run.log_scalar("timestep", t)
    _run.log_scalar("edge_density",ind.edge_density(graph))
    _run.log_scalar("strength",ind.strength(graph))
    _run.log_scalar("clustering_coeff",ind.clustering_coeff(graph))
    _run.log_scalar("nn_degree",ind.nn_degree(graph))
    _run.log_scalar("number_conflicts",len(bs.traf.cd.confpairs_all))

    comp_confs = ind.comp_conf(graph)
    _run.log_scalar("number_comp_conf",len(comp_confs))

    conf_lengths = [len(comp_conf) for comp_conf in comp_confs]
    if conf_lengths:
        conf_size = max(conf_lengths)
    else:
        conf_size = 0
        
    _run.log_scalar("conf_size", conf_size)


    

def append_variables(variables, graph):
    """
    During a time window when conflicts are present stores the values of the variables
    in the dictionary so maximum values can be computed later
    """
    variables["ed"].append(ind.edge_density(graph))
    variables["s"].append(ind.strength(graph))
    variables["cc"].append(ind.clustering_coeff(graph))
    variables["nnd"].append(ind.nn_degree(graph))

    for pair in bs.traf.cd.confpairs_all:
        variables["confs"].append({pair[0], pair[1]})

    comp_confs = ind.comp_conf(graph)
    for conf in comp_confs:
        variables["comp_confs"].append(conf)

def log_conflict_variables(t, variables, _run):
    _run.log_scalar("timestep", t)

    _run.log_scalar("edge_density", max(variables["ed"]))
    _run.log_scalar("strength",max(variables["s"]))
    _run.log_scalar("clustering_coeff",max(variables["c"]))
    _run.log_scalar("nn_degree",max(variables["nnd"]))

    _run.log_scalar("number_conflicts", len(variables["confs"]))
    _run.log_scalar("number_comp_conf",len(variables["comp_confs"]))

    conf_size = max([len(comp_conf) for comp_conf in variables["comp_confs"]])
    _run.log_scalar("conf_size", conf_size)


@ex.config
def cfg():
    
    center = (47, 9)
    radius = 0.001
    n_ac = 500
    sim_time =1*60
    n_sources = 10
    n_runs = 1

@ex.automain
def complexity_simulation(_run, center, radius, n_ac, sim_time, n_sources, n_runs):
    

    bs.init('sim-detached')
    ## We initialize the simulation ##

    for run in range(n_runs):

        bs.sim.simdt = 1
        bs.sim.simt = 0
        t_max = sim_time #15 mins

        ntraf = bs.traf.ntraf
        n_steps = t_max//bs.sim.simdt + 1
        t = np.linspace(0, t_max, n_steps)

        sources_position = create_sources(center, radius, n_sources)
        init_at(radius, center, n_ac)

        variables = {"ed": [], "s": [], "cc": [], "nnd": [], "confs": [], "comp_confs": []}
        """ Main loop """
        change_ffmode()
        for i in range(n_steps):

            
            """ Check if the acs are out of bounds and delete them if so """
            check_boundaries(bs.traf, center, radius)

            """ Spawning aircrafts in the sources """
            if bs.traf.ntraf < n_ac:
                spawn_ac(sources_position, radius, center, number_of_aircrafts = n_ac - bs.traf.ntraf)

            graph = at_to_graph()

            if (len(bs.traf.cd.confpairs_all) == 0) | (i == n_steps - 1):

                if sum([len(variables[key]) for key in variables]) != 0:
                    log_conflict_variables(bs.sim.simt, variables, _run)

                    for key in variables:
                        variables[key] = [] # reset all the values to an empty list

                log_variables(bs.sim.simt, graph, _run)

            else:
                append_variables(variables, graph)

            #plot_at(center, radius, sources_position)
            print(len(bs.traf.cd.confpairs_all))
            simstep()
            

        bs.sim.reset()
    
    
