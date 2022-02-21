import numpy as np
import bluesky as bs
from bluesky import traffic as tr
from bluesky import settings
from bluesky.traffic.route import Route
from bluesky.navdatabase import Navdatabase
from bluesky.simulation import Simulation
from bluesky.traffic.performance.perfbase import PerfBase
import matplotlib.pyplot as plt
from datalogger import Logger
import rand_conflict
from bluesky.tools import geo
import random



class ScreenDummy:

 	def __init__(self):
 		pass

 	def echo(self, text="", flags = 0):
 		pass

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

	dist_to_center = 3*radius/4
	earth_radius = geo.rwgs84(center[0])/1852

	for i in range(n_sources):
		alpha = 2*np.pi/n_sources * i
		x, y = dist_to_center * np.cos(alpha), dist_to_center * np.sin(alpha)
		sources_positions.append([center[0] + 180/np.pi * y/earth_radius, center[1] + 180/np.pi * x/earth_radius])

	return sources_positions

def spawn_acs(traf, sources_position, radius, center):

	alpha = np.linspace(0, 2*np.pi, 1000)
	y = radius * np.sin(alpha)
	x = radius * np.cos(alpha)
	earth_radius = geo.rwgs84(center[0])/1852

	lats = center[0] + 180/np.pi * y/earth_radius
	lons = center[1] + 180/np.pi * x/earth_radius

	bound_coordinates = [[lat, lon] for lat, lon in zip(lats, lons)]

	for i, source in enumerate(sources_position):

		random_point = random.choice(bound_coordinates)
		while geo.latlondist(random_point[0], random_point[1], source[0], source[1]) < radius *1852:
			random_point = random.choice(bound_coordinates)

		acid = str(random.getrandbits(32))
		
		heading = geo.qdrdist(source[0], source[1], random_point[0], random_point[1])[0]

		traf.cre(acid, actype="M200", aclat=source[0], aclon=source[1], acspd=40, achdg=heading)

		#print(random_point[0], random_point[1])

		#wpname = "wp" + acid
		#bs.stack.stack('DEST {0} {1} {2} {3}'.format(acid, wpname, random_point[0], random_point[1]))
		#print(random_point, acid)
		#Route.addwpt(iac = acid, name = acid, wptype = "DEST", lat = random_point[0], lon = random_point[1])


def plot_at(center, radius, sources_position):

	alpha = np.linspace(0, 2*np.pi, 1000)
	ybound = radius * np.sin(alpha)
	xbound = radius * np.cos(alpha)
	earth_radius = geo.rwgs84(center[0])/1852

	lats = center[0] + 180/np.pi * ybound/earth_radius
	lons = center[1] + 180/np.pi * xbound/earth_radius

	x = bs.traf.lat
	y = bs.traf.lon
	vx = bs.traf.gseast
	vy = bs.traf.gsnorth

	plt.figure(figsize=(8, 8))
	plt.scatter(center[0], center[1], 2)
	plt.scatter(lats, lons, 1)
	plt.quiver(x, y, vy, vx)
	plt.show()


def complexity_simulation(ScreenDummy, center, radius, sim_time, n_sources, ac_freq):
	# Initialize global settings
	settings.init("")
	# Manually set the performance model to the one defined in the settings before
	PerfBase.setdefault(settings.performance_model)
	# Init dummy screen
	bs.scr = ScreenDummy()
	# Manually create singletons

	traf = tr.Traffic()
	bs.traf = traf

	navdb = Navdatabase()
	bs.navdb = navdb

	sim = Simulation()
	bs.sim = sim


	## We initialize the simulation ##

	bs.sim.simdt = 1
	bs.sim.simt = 0
	t_max = sim_time #15 mins

	ntraf = traf.ntraf
	n_steps = t_max//bs.sim.simdt + 1
	t = np.linspace(0, t_max, n_steps)

	sources_position = create_sources(center, radius, n_sources)
	
	logger = Logger("lat", "lon", "alt", dt = 10, name = "COMP_LOGGER")
	
	""" Main loop """
	for i in range(n_steps):

		
		logger.log()
		
		""" Check if the acs are out of bounds and delete them if so """
		check_boundaries(traf, center, radius)

		""" Spawning aircrafts in the sources """
		if bs.sim.simt % ac_freq == 0:
			spawn_acs(traf, sources_position, radius, center)


		bs.sim.simt += bs.sim.simdt

		traf.update()
		#plot_at(center, radius, sources_position)
		

	logger.stop()
	del logger

if __name__ == '__main__':
	complexity_simulation(ScreenDummy, (47, 9), 1, 15 * 60, 100, 10)

	