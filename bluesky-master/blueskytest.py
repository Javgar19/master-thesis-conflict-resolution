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

class ScreenDummy:

 	def __init__(self):
 		pass

 	def echo(self, text="", flags = 0):
 		pass


def main(ScreenDummy):
	# Initialize global settings
	settings.init("")
	# Manually set the performance model to the one defined in the settings before
	PerfBase.setdefault(settings.performance_model)
	# Init dummy screen
	bs.scr = ScreenDummy()
	# Manually create singletons
	n = 4
	traf = rand_conflict.generate(n, 2000, 3000, t_loss = 15, t_la = 8)
	bs.traf = traf

	navdb = Navdatabase()
	bs.navdb = navdb

	sim = Simulation()
	bs.sim = sim

	## We create the planes ##
	n = 4




	## We initialize the simulation ##

	bs.sim.simdt = 1
	bs.sim.simt = 0
	t_max = 4000

	ntraf = traf.ntraf
	n_steps = t_max//bs.sim.simdt + 1
	t = np.linspace(0, t_max, n_steps)

	res = np.zeros((n_steps, 4, ntraf))
	
	logger = Logger("lat", "lon", "alt", dt = 100, name = "SECOND_LOGGER", aircrafts_id=["001"])
	
	for i in range(n_steps):
		logger.log()
		res[i] = [bs.traf.lat, bs.traf.lon, bs.traf.alt, bs.traf.tas]

		bs.sim.simt += bs.sim.simdt

		traf.update()

	logger.stop()
	del logger
	
	plt.plot(t, res[:, 3, bs.traf.id.index(acids[1])])
	plt.show()
	# print(bs.traf.id)

if __name__ == '__main__':
	main(ScreenDummy)