from multiprocessing.dummy import Process
import sys
import numpy as np

import bluesky as bs
from bluesky import traffic as tr
from bluesky import settings
from bluesky.stack.simstack import process
from bluesky.traffic.route import Route
from bluesky.navdatabase import Navdatabase
from bluesky.simulation import Simulation
from bluesky.traffic.performance.perfbase import PerfBase
import matplotlib.pyplot as plt
from bluesky.tools import geo
import random
from sacred import Experiment
import networkx as nx
from multiprocessing import Lock, Process, Queue, current_process,cpu_count
import queue
import time
import complexity_indicators as ind

def prog_bar(total, progress):
    
    barLength, status = 50, ""
    num_done = progress
    progress = float(progress) / float(total)
    if progress > 1.:
        return
    block = int(round(barLength * progress))
    text = "\r[{}] {:.0f}% {}".format(
        "#" * block + "-" * (barLength - block), round(progress * 100, 0),
        status)
    sys.stdout.write(text)
    sys.stdout.flush()
 


def worker(n_runs,task_queue, done):
    bs.init('sim-detached')
    while True:
        try:
            task = task_queue.get_nowait()
        except queue.Empty:
            
            break
        else:
            answer = task()
            done.put(answer)
            #prog_bar(n_runs, done.qsize())
            time.sleep(.5)
            
            


class Simulation(object):
    def __init__(self, center, radius, n_ac, sim_time, rpz, tcpa_thresh, num_sim) -> None:
        self.center = center
        self.radius = radius
        self.n_ac = n_ac
        self.sim_time = sim_time
        self.rpz = rpz
        self.tcpa_thresh = tcpa_thresh
        self.num_sim = num_sim

    def __call__(self):

        ntraf = bs.traf.ntraf
        bs.traf.cd.setmethod("ON")
        bs.traf.cd.rpz_def = self.rpz
        bs.traf.cd.dtlookahead_def = 15
        init_at(self.radius, self.center, self.n_ac)

        variables = {"ed": [], "s": [], "cc": [], "nnd": [], "confs": [], "comp_confs": []}
        
        """ Main loop """
        t = 0
        result = []
        while bs.sim.simt <=self.sim_time:
            
            if not bs.sim.ffmode:
                change_ffmode()


            """ Check if the acs are out of bounds and delete them if so """
            check_boundaries(bs.traf, self.center, self.radius)

            """ Spawning aircrafts """
            if bs.traf.ntraf < self.n_ac:
                spawn_ac(self.radius, self.center, number_of_aircrafts = self.n_ac - bs.traf.ntraf)

            graph = at_to_graph(4, self.tcpa_thresh)

            if (len(bs.traf.cd.confpairs_unique) == 0) | (bs.sim.simt + bs.sim.simdt >= self.sim_time):

                if sum([len(variables[key]) for key in variables]) != 0:
                    result.append(log_conflict_variables(bs.sim.simt, variables, self.num_sim))
                    return result

                    for key in variables:
                        variables[key] = [] # reset all the values to an empty list

                
                result.append(log_variables(bs.sim.simt, graph, self.num_sim))

            else:
                append_variables(variables, graph)

            simstep()
            t = bs.sim.simt
            print(self.num_sim, len(bs.traf.cd.confpairs_unique))
        
        bs.sim.reset()
        return result


    def __str__(self) -> str:
        return f"Run number {self.num_sim}"



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


def simstep():
    bs.sim.step()
    bs.net.step()
    

def change_ffmode(mode=True):
        bs.sim.ffmode = mode
        bs.sim.dtmult = 50.0


def goal_state(source_pos, radius,center):
    distance_to_goal = radius * np.sqrt(random.random())
    
    angle = random.random() * 360.
    goal_lat, goal_lon = geo.qdrpos(center[0], center[1], angle, distance_to_goal/2)
    return goal_lat, goal_lon

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
    
    if len(tcpa_matrix) != 0:
        tcpa = np.abs(tcpa_matrix[bs.traf.id.index(acid1)][bs.traf.id.index(acid2)])
        
    else:
        tcpa = 0

    if tcpa <= Htcp:
        return 1

    elif tcpa >= threshtcp:
        return 0

    else:
        return (threshtcp - tcpa)/(threshtcp) 

def at_to_graph(H, thresh):
    graph = nx.Graph()
    graph.add_nodes_from(bs.traf.id)
    for pair in bs.traf.cd.confpairs_unique:
        pair = list(pair)
        wtcp = tcp(pair[0], pair[1], H, thresh)
        graph.add_edge(pair[0], pair[1], weight = wtcp)

    return graph

def log_variables(t, graph, num_sim):
    result = []
    result.append(t)
    result.append(num_sim)
    result.append(ind.edge_density(graph))
    result.append(ind.strength(graph))
    result.append(ind.clustering_coeff(graph))
    result.append(ind.nn_degree(graph))
    result.append(len(bs.traf.cd.confpairs_unique))

    comp_confs = ind.comp_conf(graph)
    result.append(len(comp_confs))

    conf_lengths = [len(comp_conf) for comp_conf in comp_confs]
    if conf_lengths:
        conf_size = max(conf_lengths)
    else:
        conf_size = 0
        
    result.append(conf_size)
    return result



def append_variables(variables, graph):
    """
    During a time window when conflicts are present stores the values of the variables
    in the dictionary so maximum values can be computed later
    """
    variables["ed"].append(ind.edge_density(graph))
    variables["s"].append(ind.strength(graph))
    variables["cc"].append(ind.clustering_coeff(graph))
    variables["nnd"].append(ind.nn_degree(graph))

    for pair in bs.traf.cd.confpairs_unique:
        pair = list(pair)
        pair = {pair[0], pair[1]}
        if (not pair in variables["confs"]):
            variables["confs"].append(pair)

    comp_confs = ind.comp_conf(graph)
    for conf in comp_confs:
        if not conf in variables["comp_confs"]:
            variables["comp_confs"].append(conf)

def log_conflict_variables(t, variables, num_sim):
    result = []
    result.append(t)
    result.append(num_sim)

    result.append(max(variables["ed"]))
    result.append(max(variables["s"]))
    result.append(max(variables["cc"]))
    result.append(max(variables["nnd"]))

    result.append(len(variables["confs"]))
    result.append(len(variables["comp_confs"]))


    if variables["comp_confs"]:
        conf_size = max([len(comp_conf) for comp_conf in variables["comp_confs"]])
    else:
        conf_size = 0

    result.append(conf_size)
    return result

def log_after_par(results, _run):
    _run.log_scalar("timestep", results[0])
    _run.log_scalar("num_sim", results[1])
    _run.log_scalar("edge_density", results[2])
    _run.log_scalar("strength", results[3])
    _run.log_scalar("clustering_coeff", results[4])
    _run.log_scalar("nn_degree", results[5])
    _run.log_scalar("number_conflicts", results[6])
    _run.log_scalar("number_comp_conf", results[7])
    _run.log_scalar("conf_size",  results[8])


@ex.config
def cfg():
    
    center = (47, 9)
    radius = 0.5
    n_ac = 100
    sim_time = 2*radius*1850*0.1
    n_runs = 1
    rpz = 0.089
    tcpa_thresh = 35

@ex.automain
def complexity_simulation(_run, center, radius, n_ac, sim_time, n_runs, rpz, tcpa_thresh):

    runs_to_do = Queue()
    runs_done = Queue()
    procs = []
    res = []
    
    
   
    n_runs += 1 #queue.get_nowait() ends prematurely, need to figure out why
    for i in range(n_runs):
        runs_to_do.put(Simulation(center, radius, n_ac, sim_time, rpz, tcpa_thresh, i))

    num_workers = cpu_count()
    print(f"{num_workers} PROCESSES ARE BEING CREATED")
    for i in range(num_workers):
        p = Process(target=worker,args=(n_runs,runs_to_do,runs_done))
        procs.append(p)
        p.start()

    for p in procs:
        print("JOINING")
        p.join()
    
    while not runs_done.empty():
        res.append(runs_done.get())

    print(res)
    
    # We log with Incense every run after the process has finished #
    for run in res:
        log_after_par(run, _run)
    
    