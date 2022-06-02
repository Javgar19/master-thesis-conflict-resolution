import random
import math
import numpy as np
import bluesky as bs
from bluesky import traffic as tr
from bluesky.traffic.asas import statebased
from bluesky.tools import geo
from bluesky.tools.aero import vcasormach, nm, casormach2tas, tas2cas, ft
from bluesky.traffic.windfield import Windfield
from bluesky.traffic.autopilot import Autopilot
import networkx as nx
import matplotlib.pyplot as plt


class Aircraft:
	def __init__(self, acid, lat, lon, spd, hdg, alt = 0):
		self.acid = acid
		self.lat = lat
		self.lon = lon 
		self.spd = spd
		self.hdg = hdg 
		self.alt = alt

		self.trk = hdg 

		#velocities
		self.tas, self.cas, self.M = vcasormach(spd, alt)
		self.gs = self.tas 
		hdgrad = np.radians(hdg)
		self.gsnorth = self.tas * np.cos(hdgrad)
		self.gseast = self.tas * np.sin(hdgrad)
		self.vs = 0

def at_var(at: [Aircraft], var):
	return np.array([getattr(ac, var) for ac in at])

def plot_at(at: [Aircraft]):
	lats, lons = at_var(at, "lat"), at_var(at, "lon")
	vx, vy = at_var(at, "gseast"), at_var(at, "gsnorth")

	plt.quiver(lons, lats, vx, vy)
	plt.title(f"Conflict scenario with {len(at)} aircrafts")
	plt.xlabel("Longitude")
	plt.ylabel("Latitude")
	plt.show()


def generate(target, spd_min, spd_max, t_loss, t_la):
	"""
	Returns a conflict scenario between a certain number of aircrafts in the form of a Traffic object

	Arguments
	-----------------------------
		-target: number of aircraft involved in the conflict
		-spd_min, spd_max = aircraft speed boundaries
		-t_loss: time till conflict
		-t_la: look-ahead time to avoid accidental early conflicts
	"""
	
	created = []

	hdgs_conflict = [0, 45, 90, 135, 180, -45, -135, -90]
	var = 10
	lat_ref = 41.4
	lon_ref = 2.15
	spd_ref = random.uniform(spd_min, spd_max)
	hdg_ref = random.uniform(1, 360)
	
	created.append(Aircraft(acid = "ac_ref", 
				lat=lat_ref, lon=lon_ref,
				spd=spd_ref, hdg=hdg_ref))

	while (len(created) < target):
		accepted = False

		while accepted == False:
			hdg = random.choice(hdgs_conflict) + random.uniform(-var, var)
			severity = random.uniform(0.0, 1.0)
			cpa = 240 - 240 * severity

			chosen_ac = random.choice(created)

			# create the chosen ac in ac_proposed traffic object 
			proposed_acs = created.copy()
			
			# create ac in conflict with the chosen one
			ac_proposed_id = "ac_" + str(len(created) + 1)
			ac_proposed = creconfs(acid = ac_proposed_id,
							 	 target = chosen_ac, dpsi = hdg,
							 	 dcpa = cpa, tlosh = t_loss)

			# check for accidental conflicts in one look-ahead time
			proposed_acs.append(ac_proposed)
			created.append(ac_proposed)
			in_conf = detect(ownship = created, 
										 intruder = proposed_acs, 
										 rpz = 2, hpz = 3, dtlookahead = t_la)

			if not in_conf[0]:
				accepted = True
				print("+1!!!!")
				# add the conflict ac to the created trafic object
			else:
				created.drop(ac_proposed)

	return created

def creconfs(acid, target, dpsi, dcpa, tlosh, dH=None, tlosv=None, spd=None):
    ''' Create an aircraft in conflict with target aircraft.

        Arguments:
        - acid: callsign of new aircraft
        - actype: aircraft type of new aircraft
        - targetidx: id (callsign) of target aircraft
        - dpsi: Conflict angle (angle between tracks of ownship and intruder) (deg)
        - cpa: Predicted distance at closest point of approach (NM)
        - tlosh: Horizontal time to loss of separation ((hh:mm:)sec)
        - dH: Vertical distance (ft)
        - tlosv: Vertical time to loss of separation
        - spd: Speed of new aircraft (CAS/Mach, kts/-)
    '''
    latref  = target.lat  # deg
    lonref  = target.lon  # deg
    altref  = target.alt  # m
    trkref  = np.radians(target.trk)
    gsref   = target.gs   # m/s
    tasref  = target.tas   # m/s
    vsref   = target.vs   # m/s
    cpa     = dcpa * nm
    pzr     = bs.settings.asas_pzr * nm
    pzh     = bs.settings.asas_pzh * ft
    trk     = trkref + np.radians(dpsi)

    if dH is None:
        acalt = altref
        acvs  = 0.0
    else:
        acalt = altref + dH
        tlosv = tlosh if tlosv is None else tlosv
        acvs  = vsref - np.sign(dH) * (abs(dH) - pzh) / tlosv

    if spd:
        # CAS or Mach provided: convert to groundspeed, assuming that
        # wind at intruder position is similar to wind at ownship position
        tas = tasref if spd is None else casormach2tas(spd, acalt)
        tasn, tase = tas * np.cos(trk), tas * np.sin(trk)
        #wn, we = Windfield.getdata(latref, lonref, acalt)
        gsn, gse = tasn, tase 
    else:
        # Groundspeed is the same as ownship
        gsn, gse = gsref * np.cos(trk), gsref * np.sin(trk)

    # Horizontal relative velocity vector
    vreln, vrele = gsref * np.cos(trkref) - gsn, gsref * np.sin(trkref) - gse
    # Relative velocity magnitude
    vrel    = np.sqrt(vreln * vreln + vrele * vrele)
    # Relative travel distance to closest point of approach
    drelcpa = tlosh * vrel + (0 if cpa > pzr else np.sqrt(pzr * pzr - cpa * cpa))
    # Initial intruder distance
    dist    = np.sqrt(drelcpa * drelcpa + cpa * cpa)
    # Rotation matrix diagonal and cross elements for distance vector
    rd      = drelcpa / dist
    rx      = cpa / dist
    # Rotate relative velocity vector to obtain intruder bearing
    brn     = np.degrees(math.atan2(-rx * vreln + rd * vrele,
                             rd * vreln + rx * vrele))

    # Calculate intruder lat/lon
    aclat, aclon = geo.kwikpos(latref, lonref, brn, dist / nm)
    # convert groundspeed to CAS, and track to heading using actual
    # intruder position
    #wn, we     = Windfield.getdata(aclat, aclon, acalt)
    tasn, tase = gsn, gse 
    acspd      = tas2cas(np.sqrt(tasn * tasn + tase * tase), acalt)
    achdg      = np.degrees(math.atan2(tase, tasn))

    # Create and, when necessary, set vertical speed
   #Autopilot.selaltcmd(len(self.lat) - 1, altref, acvs)
   #self.vs[-1] = acvs
    return Aircraft(acid, aclat, aclon, acspd, achdg, acalt)



def detect(ownship, intruder, rpz, hpz, dtlookahead):
    ''' Conflict detection between ownship (traf) and intruder (traf/adsb).'''
    # Identity matrix of order ntraf: avoid ownship-ownship detected conflicts
    I = np.eye(len(ownship))

    # Horizontal conflict ------------------------------------------------------

    # qdrlst is for [i,j] qdr from i to j, from perception of ADSB and own coordinates
    qdr, dist = geo.kwikqdrdist_matrix(np.asmatrix(at_var(ownship,"lat")), np.asmatrix(at_var(ownship,"lon")),
                                np.asmatrix(at_var(intruder,"lat")), np.asmatrix(at_var(intruder,"lon")))

    # Convert back to array to allow element-wise array multiplications later on
    # Convert to meters and add large value to own/own pairs
    qdr = np.asarray(qdr)
    dist = np.asarray(dist) * nm + 1e9 * I

    # Calculate horizontal closest point of approach (CPA)
    qdrrad = np.radians(qdr)
    dx = dist * np.sin(qdrrad)  # is pos j rel to i
    dy = dist * np.cos(qdrrad)  # is pos j rel to i

    # Ownship track angle and speed
    owntrkrad = np.radians(at_var(ownship,"trk"))
    ownu = at_var(ownship,"gs") * np.sin(owntrkrad).reshape((1, len(ownship)))  # m/s
    ownv = at_var(ownship,"gs") * np.cos(owntrkrad).reshape((1, len(ownship)))  # m/s

    # Intruder track angle and speed
    inttrkrad = np.radians(at_var(intruder,"trk"))
    intu = at_var(intruder,"gs") * np.sin(inttrkrad).reshape((1, len(ownship)))  # m/s
    intv = at_var(intruder,"gs") * np.cos(inttrkrad).reshape((1, len(ownship)))  # m/s

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
    dalt = at_var(ownship,"alt").reshape((1, len(ownship))) - \
        at_var(intruder,"alt").reshape((1, len(ownship))).T  + 1e9 * I

    dvs = at_var(ownship, "vs").reshape(1, len(ownship)) - \
        at_var(intruder,"vs").reshape(1, len(ownship)).T
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
    tcpamax = np.max(tcpa * swconfl, 1)

    # Select conflicting pairs: each a/c gets their own record
    confpairs = [(ownship[i].id, ownship[j].id) for i, j in zip(*np.where(swconfl))]
    swlos = (dist < rpz) * (np.abs(dalt) < hpz)
    lospairs = [(ownship[i].id, ownship[j].id) for i, j in zip(*np.where(swlos))]

    return confpairs, lospairs, inconf, tcpamax, \
        qdr[swconfl], dist[swconfl], np.sqrt(dcpa2[swconfl]), \
            tcpa[swconfl], tinconf[swconfl]




#	def generate(target, spd_min, spd_max, t_loss, t_la):
# 	created = list()
# 	hdgs_conflict = [0, 45, 90, 135, 180, -45, -135, -90]
# 	var = 10
# 	lat_ref = 41.4
# 	lon_ref = 2.15
# 	spd_ref = random.uniform(spd_min, spd_max)
# 	hdg_ref = random.uniform(1, 360)

# 	bs.traf.cre(acid = "ac_ref" + str(len(created)), # no ac object wtf
# 						 aclat=lat_ref, aclon=lon_ref,
# 						 acspd=spd_ref, achdg=hdg_ref)

# 	while (len(bs.traf.id) < target):
# 		accepted = False

# 		while accepted == False:
# 			hdg = random.sample(hdgs_conflict, 1) + random.uniform(-var, var)
# 			severity = random.uniform(0.0, 1.0)
# 			cpa = 240 - 240 * severity
# 			chosen = random.sample(bs.traf.id, 1)

# 			ac_proposed_id = "ac_" + str(len(bs.traf.id) + 1)
# 			bs.traf.creconfs(acid = ac_proposed_id, actype = "B744",
# 							 targetidx = bs.traf.id.index(chosen), dpsi = hdg,
# 							 dcpa = cpa, tlosh = t_loss)
# 			statebased.StateBased.detect(ownship = bs.traf, 
# 										 intruder = ac_proposed_id, 
# 										 rpz, hpz, dtlookahead = t_la)