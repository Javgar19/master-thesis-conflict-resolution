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
import networkx as nx



def simstep():
    bs.sim.step()
    bs.net.step()

def change_ffmode(mode=True):
    bs.sim.ffmode = mode
    bs.sim.dtmult = 5.0

def plot_at():
    x = bs.traf.lat
    y = bs.traf.lon
    vx = bs.traf.gseast
    vy = bs.traf.gsnorth

    plt.figure(figsize=(8, 8))
    plt.quiver(x, y, vy, vx)
    plt.title("Scenario initial configuration")
    plt.show()

def dhij(acid1, acid2):
    """
    Horizontal distance (m) between acid1 and acid2
    """
    return geo.latlondist(bs.traf.lat[bs.traf.id.index(acid1)], bs.traf.lon[bs.traf.id.index(acid1)],
                        bs.traf.lat[bs.traf.id.index(acid2)], bs.traf.lon[bs.traf.id.index(acid2)])

def whij(acid1, acid2, H, thresh):
    """
    Horizontal weight for the edge between acid1 and acid2 
    """
    dij = dhij(acid1, acid2)

    if dij <= H:
        return 1

    elif dij >= thresh:
        return 0

    else:
        return (thresh - dij)/(thresh)

def tcp(acid1, acid2, Htcp, threshtcp):
    """
    Horizontal weight for the edge between acid1 and acid2 
    """
    if len(bs.traf.cd.tcpa) != 0:
        tcpa = np.abs(bs.traf.cd.tcpa[0][1])
    else:
        tcpa = 0

    if tcpa <= Htcp:
        return 1

    elif tcpa >= threshtcp:
        return 0

    else:
        return (threshtcp - tcpa)/(threshtcp)

def tcpa_simulation():
    
    bs.init('sim-detached')
    ## We initialize the simulation ##

    bs.traf.cd.setmethod("ON")
    bs.traf.cd.rpz_def = 5
    bs.traf.cd.dtlookahead_def = 300
    #bs.traf.cd.rpz = rpz
    #bs.traf.cd.dtlookahead = 300

 

    #bs.sim.simdt = 1
    #bs.sim.simt = 0
    t_max = 60 #15 mins

    ntraf = bs.traf.ntraf
    n_steps = int(t_max//bs.sim.simdt + 1)
    t = np.linspace(0, t_max, n_steps)

    position1 = (47.0, 9.0)
    angle = 0
    distance = 200/1852
    position2 = geo.qdrpos(position1[0], position1[1], angle, distance)
    velocity = 40

    H = 0
    thresh = 200
    Htcp = 0
    threshtcp = 2.50167139e+00

    bs.traf.cre("001", actype="M200", aclat=position1[0], aclon=position1[1], acspd=velocity, achdg=angle)
    bs.traf.cre("002", actype="M200", aclat=position2[0], aclon=position2[1], acspd=velocity, achdg=angle + 180)
    """ Main loop """
    change_ffmode()

    wdistance = []
    wtcpa = []

    for i in range(n_steps):
        wdistance.append(whij("001", "002", H, thresh))
        wtcpa.append(tcp("001", "002", Htcp, threshtcp))

        
      
        
        simstep()
        

    bs.sim.reset()
    plt.plot(wtcpa[1:150], label = "tcpa")
    plt.plot(wdistance[0:150], label = "distance")
    plt.legend()
    plt.xlabel("Iteration")
    plt.ylabel("Weight")
    plt.show()

    print(wtcpa[10:50:10])
    print(wdistance[10:50:10])

if __name__=='__main__':
    tcpa_simulation()

