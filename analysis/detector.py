##################################################################################################
##################################################################################################

#THIS FILE CONTAINS DETECTOR GEOMETRY INFORMATION#

##################################################################################################

import numpy as np

class Detector:

	BoxLimits = [  [-5000.0, 5000.0],  [6000.0, 8917.0],  [7000.0, 17000.0]    ]
	LayerYLims = [ [6001., 6004.],  [6104., 6107.], [8001., 8004.], [8104., 8107.], [8501., 8504.], [8604., 8607.], [8707., 8710.], [8810., 8813.], [8913., 8916.]  ]
	ModuleXLims = [ [-4950. + 1000.*n, -4050. + 1000*n] for n in range(10) ]
	ModuleZLims = [ [7000.  + 1000.*n,  7900. + 1000*n] for n in range(10) ]

	#               x range              y range             z range

	def __init__(self):
		# print("Detector Constructed")
		pass
	def xLims(self):
		return self.BoxLimits[0]

	def yLims(self):
		return self.BoxLimits[1]

	def zLims(self):
		return self.BoxLimits[2]

	def LayerY(self, n):
		return self.LayerYLims[n]

	def LayerYMid(self, n):
		return (self.LayerYLims[n][0] + self.LayerYLims[n][1])/2.

	def numLayers(self):
		return len(self.LayerYLims)

	def DrawColor(self):
		return "tab:gray"

	def inLayer(self, yVal):
		for layerN, layerLims in enumerate(self.LayerYLims):
			if yVal > layerLims[0] and yVal < layerLims[1]:
				return layerN
		return -1

	def nextLayer(self, yVal):
		for n in range(len(self.LayerYLims)-1):
			if yVal > self.LayerYLims[n][1] and yVal < self.LayerYLims[n+1][0]:
				return n+1
		return 999


	def inLayer_w_Error(self, yVal, yErr):
		for layerN, layerLims in enumerate(self.LayerYLims):

			lower = yVal - yErr
			upper = yVal + yErr

			if lower < layerLims[0] and upper > layerLims[1]:
				return layerN

			if lower < layerLims[0] and upper > layerLims[0]:
				return layerN

			if lower < layerLims[1] and upper > layerLims[1]:
				return layerN


		return -1

	def inBox(self, x, y, z):
		if x > self.xLims()[0] and x < self.xLims()[1]:
			if y > self.yLims()[0] and y < self.yLims()[1]:
				if z > self.zLims()[0] and z < self.zLims()[1]:
					return True
		return False


	#determine number of layers a track goes through
	def nLayers(self, x0, y0, z0, vx, vy, vz):
		count = 0
		for n in range(len(self.LayerYLims)):
			layerY = self.LayerYMid(n)
			if (layerY-y0)/vy < 0:
				continue
			else:
				dt = (layerY - y0)/vy

				x1 = x0 + dt*vx
				z1 = y0 + dt*vz

				if inBox(x1, layerY, z1):
					count += 1

		return count

	#determine number of SENSITIVE layers a track goes through
	def nSensitiveLayers(self, x0, y0, z0, vx, vy, vz):
		count = 0
		for n in range(len(self.LayerYLims)):
			layerY = self.LayerYMid(n)
			if (layerY-y0)/vy < 0:
				continue
			else:
				dt = (layerY - y0)/vy

				x1 = x0 + dt*vx
				z1 = y0 + dt*vz

				if self.inSensitiveElement(x1, layerY, z1):
					count += 1

		return count


	##get points inside the detector for reconstructed track
	def RecoTrackPoints(self, x0, y0, z0, vx, vy, vz):
		x, y, z = [], [], []
		_x, _y, _z = x0, y0, z0

		time_spacing = 0.1 #ns

		while self.inBox(_x, _y, _z):
			x.append(_x)
			y.append(_y)
			z.append(_z)

			_x += vx*time_spacing
			_y += vy*time_spacing
			_z += vz*time_spacing


		return x, y, z


	def FindIntercept(self, x0, y0, z0, vx, vy, vz):
		#x, y, z = [], [], []
		#_x, _y, _z = x0, y0, z0


		pos = np.array([x0, y0, z0])
		v = np.array([vx, vy, vz])
		t = []

		#print("v is {}".format(v))

		for i in range(len(pos)):
			if v[i] > 0:
				t.append((self.BoxLimits[i][1] - pos[i]) / v[i])

			else:
				t.append((self.BoxLimits[i][0] - pos[i]) / v[i])

		t_intercept = np.amin(t)
		#print("the ts are {}".format(t))
		#print("t intercept is {}".format(t_intercept))

		return (pos + t_intercept * v)


	def inModuleX(self, xVal):
		for moduleN, moduleLims in enumerate(self.ModuleXLims):
			if xVal > moduleLims[0] and xVal < moduleLims[1]:
				return moduleN
		return -1

	def inModuleZ(self, zVal):
		for moduleN, moduleLims in enumerate(self.ModuleZLims):
			if zVal > moduleLims[0] and zVal < moduleLims[1]:
				return moduleN
		return -1

	def inSensitiveElement(self, x, y, z):
		if self.inLayer(y) >= 0:
			if self.inModuleX(x) >= 0:
				if self.inModuleZ(z) >= 0:
					return True
		return False
