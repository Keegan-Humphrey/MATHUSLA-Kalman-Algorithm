import numpy as np 
import ROOT as root 

c = 29.98


def colorMap(pdg):
	if (pdg == 13):
		return "c"
	if (pdg == -13):
		return "c" 
	if (pdg == 11):
		return "g"
	if (pdg == -11):
		return "m"
	return "y"   



class Vector:

	def __init__(self, x, y, z):
		self.x, self.y, self.z = x, y, z

	def magnitude(self):
		return np.sqrt( self.x**2  + self.y**2 + self.z**2   )

	def M2(self):
		return ( self.x**2  + self.y**2 + self.z**2 )

	def __add__(self, other):
		nx = self.x + other.x
		ny = self.y + other.y
		nz = self.z + other.z
		return Vector(nx, ny, nz)

	def __sub__(self, other):
		nx = self.x - other.x
		ny = self.y - other.y
		nz = self.z - other.z
		return Vector(nx, ny, nz)

	def __pow__(self, c): #SCALER MULTIPLICATION
		return Vector(self.x*c, self.y*c, self.z*c)

	def __str__(self):
		return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"

	def __repr__(self):
		return self.__str__()

class Particle:

	def __init__(self, trackID, pdgID, parentTrackID=-1):
		self.trackID = trackID
		self.pdgID = pdgID
		self.parentTrackID = parentTrackID

	def __str__(self):
		_str = "PDG: " + str(self.pdgID) +  " TrackID " + str(self.trackID) + " parentTrackID: " + str(self.parentTrackID)
		return _str

	def __repr__(self):
		return "pdg=" + str(self.pdgID) + "-TrackID=" + str(self.trackID)

	def __eq__(self, other):
		return (self.trackID == other.trackID)

	def __hash__(self):
		return self.trackID

	def color(self):
		if (self.trackID == 1):
			return "r"
		return colorMap(self.pdgID)

class TrackPoint:

	def __init__(self, time, location, energy, momentum):
		self.time = time 
		self.location = location
		self.energy = energy
		self.momentum = momentum

		m = np.sqrt(energy**2 - momentum.M2())
		gamma = energy/m
		self.velocity = (momentum ** (1.0/(m*gamma))) ** c

	def __lt__(self, other):
		return self.time < other.time
	def __le__(self, other):
		return self.time <= other.time
	def __gt__(self, other):
		return self.time > other.time
	def __ge__(self, other):
		return self.time >= other.time

	

class Track:


	def __init__(self, particle):
		self.particle = particle
		self.pointList = []

	def AddPoint(self, point):
		self.pointList.append(point)

	def PointAtTime(self, t): #return position of the track at time t
		if t < self.pointList[0].time:
			return Vector(0, 0, 0)

		if t > self.pointList[len(self.pointList)-1].time:
			return self.pointList[len(self.pointList)-1].location

		for pn1 in range(len(self.pointList)-1):
			pn2 = pn1+1
			point1 = self.pointList[pn1]
			point2 = self.pointList[pn2]

			if point1.time < t and t < point2.time:
				return point1.velocity ** (t - point1.time) + point1.location 

	def __repr__(self):
		_str = str(self.particle) + "\n"
		_str += "Energy: " + str(round(self.pointList[0].energy, 2)) + " MeV\n"
		#for point in (self.pointList):
		#	_str += "time: " + str(point.time)
		#	_str += " position: " + str(point.location)
		#	_str += "\n"
		return _str

	def Name(self, pdg):
		if pdg == 13:
			return "mu-"
		if pdg == -13:
			return "mu+"
		if pdg == 11:
			return "e-"
		if pdg == -11:
			return "e+"
		return "pdg: " + str(pdg)# + ", " + + str(round(self.pointList[0].energy, 2)) + "MeV"

	def LabelString(self):
		return self.Name(self.particle.pdgID)

	def color(self):
		return self.particle.color()

	def Range(self):
		return (self.pointList[0].location - self.pointList[len(self.pointList)-1].location).magnitude()
	
	def TimeRange(self):
		return np.absolute(self.pointList[0].time - self.pointList[len(self.pointList)-1].time)


	def __lt__(self, other):
		return other.particle.trackID < self.particle.trackID

	def v0(self):
		return self.pointList[0].velocity

	def Energy(self):
		print(round(self.pointList[0].energy, 2))
		
class RecoTrackPt:
	def __init__(self, x, y, z, c=0):
		self.x = x
		self.y = y
		self.z = z
		self.c = c

	def __repr__(self):
		_str = "location: " + "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
		return _str
