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

	dist_to_center = radius

	for i in range(n_sources):
		alpha = 360/n_sources * i
		lat, lon = geo.qdrpos(center[0], center[1], alpha, dist_to_center)
		sources_positions.append([lat, lon])

	return sources_positions

def init_at(radius, center, n_ac):

	earth_radius = geo.rwgs84(center[0])/1852

	for _ in range(n_ac):

		random_angle = random.random() * 360
		random_distance = random.random()

		random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)


		angle = geo.qdrdist(random_lat, random_lon, center[0], center[1])[0]
		limit_angle = 60.0 # 60deg when the sources are in the np

		acid = str(random.getrandbits(32))
		heading = random.uniform(angle - limit_angle, angle + limit_angle)

		bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon, acspd=40, achdg=heading)


		#wpname = "wp" + acid
		#bs.stack.stack('DEST {0} {1} {2} {3}'.format(acid, wpname, random_point[0], random_point[1]))
		#print(random_point, acid)
		#Route.addwpt(iac = acid, name = acid, wptype = "DEST", lat = random_point[0], lon = random_point[1])

def spawn_ac(sources_position, radius, center, number_of_aircrafts):

	for _ in range(number_of_aircrafts):
		source = random.choice(sources_position)
		angle = geo.qdrdist(source[0], source[1], center[0], center[1])[0]
		limit_angle = 60.0 * 180 / np.pi

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


def complexity_simulation(ScreenDummy, center, radius, n_ac, sim_time, n_sources):
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
	init_at(radius, center, n_ac)
	plot_at(center, radius, sources_position)
	
	logger = Logger("lat", "lon", "alt", dt = 10, name = "COMP_LOGGER")
	
	""" Main loop """
	for i in range(n_steps):

		
		logger.log()
		
		""" Check if the acs are out of bounds and delete them if so """
		check_boundaries(traf, center, radius)

		""" Spawning aircrafts in the sources """
		if bs.traf.ntraf < n_ac:
			spawn_ac(sources_position, radius, center, number_of_aircrafts = n_ac - bs.traf.ntraf)


		bs.sim.simt += bs.sim.simdt

		traf.update()
		
		if bs.traf.ntraf != 100:
			print(bs.traf.ntraf)

	logger.stop()
	del logger

if __name__ == '__main__':
	complexity_simulation(ScreenDummy, center=(47, 9), radius=1, n_ac=100, sim_time=15 * 60, n_sources=10)

	