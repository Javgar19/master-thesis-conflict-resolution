import csv
import os
import numpy as np 
import datetime as datetime
import bluesky as bs
from bluesky.core import varexplorer as ve
from bluesky.traffic.asas import detection
import networkx as nx
from bluesky import traffic
from bluesky.tools import geo
from bluesky.tools.aero import vcasormach, nm, casormach2tas, tas2cas, ft

class Logger():
	"""
	Data logger class

	Methods:
	-__create_fname(): Generate the name of the logged file
	--------

	"""
	def __init__(self, *variables, dt: float, name: str = "", aircrafts_id: str = None):
	
		self.dt = dt
		self.variables = [variable for variable in variables]
		self.name = name
		self.tlog = 0

		if aircrafts_id:
			self.aircrafts_id = aircrafts_id  
			self.id_dynamic = False

		else:
			self.aircrafts_id = bs.traf.id
			self.id_dynamic = True

		# checks if any acid or variable does not exists
		for aircraft_id in self.aircrafts_id:

			if aircraft_id not in bs.traf.id:
				raise ValueError(f"The aircraft {aircraft_id} does not exist") 

		self.__create_directory()
		
		fname = self.__create_fname()
		self.__open_file(fname)
		

	def __update_acids(self):
		self.aircrafts_id = bs.traf.id


	def __create_fname(self) -> str:
		"""
		Generate the name of the logged file with the format: "name_yyymmdd_hh-mm-ss.csv" 
		where name referes to the constructor argument "name"
		"""

		timestamp = datetime.datetime.now().strftime('%Y%m%d_%H-%M-%S')
		fname = self.name + "_" + timestamp + ".csv"

		return fname

	def __open_file(self, fname):
		"""
		Create and open the output file with the created name and create the csv_writer object
		"""
		self.file = open(self.directory_path + "/" + fname, "w", newline = "")

		col_header = "# aircraft_id "
		for variable in self.variables:
			col_header += variable + " "

		col_header += "conflict_pairs compound_conflicts"
		col_header += "\n"

		self.file.write(col_header)

		self.csv_writer = csv.writer(self.file)

	def __close_file(self):
		"""
		Close the output file
		"""
		self.file.close()

	def __create_directory(self):
		"""
		Create a folder to store the output file and returns its path
		"""
		if self.name:
			self.directory_path = "./" + self.name + "_output"

		else:
			self.directory_path = "./" + "output"

		if(not os.path.isdir(self.directory_path)):
			os.mkdir(self.directory_path)

	def detect(self, ownship, intruder, rpz, hpz, dtlookahead):
		''' Conflict detection between ownship (traf) and intruder (traf/adsb).'''
		# Identity matrix of order ntraf: avoid ownship-ownship detected conflicts
		I = np.eye(ownship.ntraf)

        # Horizontal conflict ------------------------------------------------------

        # qdrlst is for [i,j] qdr from i to j, from perception of ADSB and own coordinates
		qdr, dist = geo.kwikqdrdist_matrix(np.asmatrix(ownship.lat), np.asmatrix(ownship.lon),
		                       	np.asmatrix(intruder.lat), np.asmatrix(intruder.lon))

		# Convert back to array to allow element-wise array multiplications later on
		# Convert to meters and add large value to own/own pairs
		qdr = np.asarray(qdr)
		dist = np.asarray(dist) * nm + 1e9 * I

		# Calculate horizontal closest point of approach (CPA)
		qdrrad = np.radians(qdr)
		dx = dist * np.sin(qdrrad)  # is pos j rel to i
		dy = dist * np.cos(qdrrad)  # is pos j rel to i

		# Ownship track angle and speed
		owntrkrad = np.radians(ownship.trk)
		ownu = ownship.gs * np.sin(owntrkrad).reshape((1, ownship.ntraf))  # m/s
		ownv = ownship.gs * np.cos(owntrkrad).reshape((1, ownship.ntraf))  # m/s

		# Intruder track angle and speed
		inttrkrad = np.radians(intruder.trk)
		intu = intruder.gs * np.sin(inttrkrad).reshape((1, ownship.ntraf))  # m/s
		intv = intruder.gs * np.cos(inttrkrad).reshape((1, ownship.ntraf))  # m/s

		du = ownu - intu.T  # Speed du[i,j] is perceived eastern speed of i to j
		dv = ownv - intv.T  # Speed dv[i,j] is perceived northern speed of i to j

		dv2 = du * du + dv * dv
		dv2 = np.where(np.abs(dv2) < 1e-6, 1e-6, dv2)  # limit lower absolute value
		vrel = np.sqrt(dv2)

		tcpa = -(du * dx + dv * dy) / dv2 + 1e9 * I

		# Calculate distance^2 at CPA (minimum distance^2)
		dcpa2 = np.abs(dist * dist - tcpa * tcpa * dv2)

		# Check for horizontal conflict
		# RPZ can differ per aircraft, get the largest value per aircraft pair
		rpz = np.asarray(np.maximum(np.asmatrix(rpz), np.asmatrix(rpz).transpose()))
		R2 = rpz * rpz
		swhorconf = dcpa2 < R2  # conflict or not

		# Calculate times of entering and leaving horizontal conflict
		dxinhor = np.sqrt(np.maximum(0., R2 - dcpa2))  # half the distance travelled inzide zone
		dtinhor = dxinhor / vrel

		tinhor = np.where(swhorconf, tcpa - dtinhor, 1e8)  # Set very large if no conf
		touthor = np.where(swhorconf, tcpa + dtinhor, -1e8)  # set very large if no conf

		# Vertical conflict --------------------------------------------------------

		# Vertical crossing of disk (-dh,+dh)
		dalt = ownship.alt.reshape((1, ownship.ntraf)) - \
		intruder.alt.reshape((1, ownship.ntraf)).T  + 1e9 * I

		dvs = ownship.vs.reshape(1, ownship.ntraf) - \
		intruder.vs.reshape(1, ownship.ntraf).T
		dvs = np.where(np.abs(dvs) < 1e-6, 1e-6, dvs)  # prevent division by zero

		# Check for passing through each others zone
		# hPZ can differ per aircraft, get the largest value per aircraft pair
		hpz = np.asarray(np.maximum(np.asmatrix(hpz), np.asmatrix(hpz).transpose()))
		tcrosshi = (dalt + hpz) / -dvs
		tcrosslo = (dalt - hpz) / -dvs
		tinver = np.minimum(tcrosshi, tcrosslo)
		toutver = np.maximum(tcrosshi, tcrosslo)

		# Combine vertical and horizontal conflict----------------------------------
		tinconf = np.maximum(tinver, tinhor)
		toutconf = np.minimum(toutver, touthor)

		swconfl = np.array(swhorconf * (tinconf <= toutconf) * (toutconf > 0.0) *
		               np.asarray(tinconf < np.asmatrix(dtlookahead).T) * (1.0 - I), dtype=np.bool)

		# --------------------------------------------------------------------------
		# Update conflict lists
		# --------------------------------------------------------------------------
		# Ownship conflict flag and max tCPA
		inconf = np.any(swconfl, 1)

		try:
			tcpamax = np.max(tcpa * swconfl, 1)
		except ValueError:
			tcpamax = 0

		# Select conflicting pairs: each a/c gets their own record
		confpairs = [(ownship.id[i], ownship.id[j]) for i, j in zip(*np.where(swconfl))]
		swlos = (dist < rpz) * (np.abs(dalt) < hpz)
		lospairs = [(ownship.id[i], ownship.id[j]) for i, j in zip(*np.where(swlos))]

		return confpairs, lospairs, inconf, tcpamax, \
		qdr[swconfl], dist[swconfl], np.sqrt(dcpa2[swconfl]), \
		    tcpa[swconfl], tinconf[swconfl]


	def __comp_conf(self, conf_pairs):
		"""
		Return a list with all compound conficts
		"""
		g = nx.Graph()
		g.add_nodes_from(self.aircrafts_id)

		for pair in conf_pairs:
			g.add_edge(pair[0], pair[1])

		if conf_pairs:
			d = list(nx.connected_components(g)) 
			return d

		return []

	def __extract_data(self):
		"""
		Returns the current data generated by the simulation
		for the specified variables and aircrafts
		"""
		# update aircraft ids
		if self.id_dynamic:
			self.__update_acids()

		aircraft_data = []
		for variable in self.variables:

			try:
				aircraft_data.append(getattr(bs.traf, variable))

			except AttributeError:
				print(f"The variable {variable} does not exist")
				exit()

		# Extract conflict pairs and add them to each record
		conf_pairs = self.detect(ownship = bs.traf, intruder = bs.traf, rpz = 0.12959, hpz = 1, dtlookahead = 15)[0]
		
		conf_data = [len(conf_pairs)/2 for _ in range(len(self.aircrafts_id))]
		
		
			
		aircraft_data.append(conf_data)

		# Compute compound conflict and add them to the records

		comp_confs = self.__comp_conf(conf_pairs)
		comp_confs = [conf for conf in comp_confs if len(conf) > 1]
		
		comp_data = [len(comp_confs) for _ in range(len(self.aircrafts_id))]


		aircraft_data.append(comp_data)

		return aircraft_data

	def log(self):
		"""
		Writes data in the csv file
		"""
		if self.file and bs.sim.simt >= self.tlog:
			# Extract the data 
			self.__extract_data()

			# Increment the tlog for the next iteration
			self.tlog += self.dt

			# Create the new row to log and log it
			new_data = self.__extract_data()

			for aircraft_id in self.aircrafts_id:

				aircraft_index = self.aircrafts_id.index(aircraft_id)

				new_row = [bs.sim.simt, aircraft_id]
				[new_row.append(new_data[i][aircraft_index]) for i, val in enumerate(new_data)]

				self.csv_writer.writerow(new_row)
	
	def stop(self):
		"""
		Close the file and reset some atributtes
		"""
		self.__close_file()
		self.tlog = 0


