import numpy as np
import physics
import visualization
from detector import Detector
import util
import ROOT as root
from joblib import dump, load
import matplotlib.pyplot as plt

#########################################################################################
### DEFINITION OF EVENT CLASS ###########################################################


class Event:
	#visEngine = visualization.Visualizer()
#	visEngine_k = visualization.Visualizer()
	###########################################################################################################################################################################
	TRUTH_PARTICLE_E_THRESHOLD = 70 #MeV

	t0 = 0.0
	tm = 0.0

	globalTime = 0.0
	timeResolution = 1.0 #ns
	timeRangeCut = 0.1 #ns

	tracksAtGlobalTime = []

#	truthTrackList = []
	recoTrackList = []
	recoVertexList = []
	digiHitList = []
	simHitList = []

	###########################################################################################################################################################################

	###########################################################################################################################################################################
	def __init__(self, Tree, EventNumber):
		self.EventNumber = EventNumber
		self.tfile = root.TFile.Open(Tree)
		self.Tree_nm = Tree
		self.Tree = self.tfile.Get("integral_tree")
#		self.Tree = self.tfile.Get("box_run")
		self.truthTrackList = []
		#self.Tree = Tree
#		self.Tree.SetBranchStatus("*", 0)

		self.writeDirectory = './'

		self.kalman = True
		self.unused = False
		self.used = False
		self.oldvertex = False
		self.vertex = False
		self.merged = False
		self.vert_vel = True

	###########################################################################################################################################################################

	###########################################################################################################################################################################
	def ExtractTruthPhysics(self):
#		self.Tree.SetBranchStatus("*", 0)
		self.Tree.SetBranchStatus("Hit_x", 1)
		self.Tree.SetBranchStatus("Hit_y", 1)
		self.Tree.SetBranchStatus("Hit_z", 1)
		self.Tree.SetBranchStatus("Hit_particlePdgId", 1)
		self.Tree.SetBranchStatus("Hit_G4ParentTrackId", 1)
		self.Tree.SetBranchStatus("Hit_G4TrackId", 1)
		self.Tree.SetBranchStatus("Hit_particlePx", 1)
		self.Tree.SetBranchStatus("Hit_particlePy", 1)
		self.Tree.SetBranchStatus("Hit_particlePz", 1)
		self.Tree.SetBranchStatus("Hit_particleEnergy", 1)
		self.Tree.SetBranchStatus("Track_numHits")
		self.Tree.SetBranchStatus("Digi_numHits")

		self.Tree.GetEntry(self.EventNumber)


		particleSet = set()

		for hitn in range(int(self.Tree.NumHits)):
			if self.Tree.Hit_particleEnergy[hitn] > self.TRUTH_PARTICLE_E_THRESHOLD:
				particleSet.add( physics.Particle( int(self.Tree.Hit_G4TrackId[hitn]), int(self.Tree.Hit_particlePdgId[hitn]), int(self.Tree.Hit_G4ParentTrackId[hitn]) ))


		excluded_pdgs = set([])
		#print(particleSet)
		for particle in particleSet:
#			print(particle.pdgID)
			if not particle.pdgID in excluded_pdgs:
        			currentTrack = physics.Track(particle)
        			for hitn in range(int(self.Tree.NumHits)):
        				if ( int(self.Tree.Hit_G4TrackId[hitn]) == particle.trackID):
        					time = self.Tree.Hit_time[hitn]
        					location = physics.Vector(self.Tree.Hit_x[hitn], self.Tree.Hit_y[hitn], self.Tree.Hit_z[hitn])
        					energy = self.Tree.Hit_particleEnergy[hitn]
        					momentum = physics.Vector( self.Tree.Hit_particlePx[hitn], self.Tree.Hit_particlePy[hitn], self.Tree.Hit_particlePz[hitn] )
        					point = physics.TrackPoint(time, location, energy, momentum)
        					currentTrack.AddPoint(point)

			if (len(currentTrack.pointList)) <= 2:
				continue

			if currentTrack.TimeRange() < self.timeRangeCut:
				continue

			currentTrack.pointList.sort()

			self.truthTrackList.append(currentTrack)

		self.ResetTracks()
		self.truthTrackList.sort()
#		self.t0 = min([track.pointList[0].time for track in self.truthTrackList])
#		self.tm = max([track.pointList[len(track.pointList)-1].time for track in self.truthTrackList])
	###########################################################################################################################################################################


	###########################################################################################################################################################################
	###########################################################################################################################################################################


	##########################################################################################################################################################################
	def Print(self):

		visEngine_local = visualization.Visualizer()

		visEngine_local.writeDirectory = self.writeDirectory

#		print(self.tm)
#		self.TruthAtTime(self.tm)

		num_truth = 0
#		list_of_truth_hits = []
		for track in self.truthTrackList:
			x = []
			y = []
			z = []
			# print(track)
			for point in track.pointList:
				x.append(point.location.x)
				y.append(point.location.y)
				z.append(point.location.z)
#				temp.append([x, y, z] )

#			list_of_truth_hits.append(temp)
#			print(y)
			visEngine_local.TrackDisplayPoints(x, y, z, track.color(), track.LabelString())

			num_truth += 1

		# print("truth tracks: ",num_truth)
		self.num_truth = num_truth
		# if num_truth > 1: ###### remove when done!!!!!
		visEngine_local.Draw(self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'truth'))




	def ResetTracks(self):
		self.tracksAtGlobalTime = [ [] for x in range(len(self.truthTrackList))  ]
		self.globalTime = self.t0

	def StepParticles(self):
		self.globalTime += self.timeResolution

		for n, track in enumerate(self.truthTrackList):
			self.tracksAtGlobalTime[n].append(track.PointAtTime(self.globalTime))


	#returns location of particles at a given time t, called by visulaization engine
	def TruthAtTime(self, t):
		if ( np.absolute(self.globalTime - t) < self.timeResolution ):
			return self.tracksAtGlobalTime

		if (t < self.globalTime):
			self.globalTime = self.t0
			self.ResetTracks()

		while (self.globalTime < t):
			self.StepParticles()

		return self.tracksAtGlobalTime

	#draw reconstructed tracks in detector det
	def DrawReco(self):
		list_of_trackPt_lists = []
		list_of_colors = []
		det = Detector()

		self.Tree.SetBranchStatus("Track_x0", 1)
		self.Tree.SetBranchStatus("Track_y0", 1)
		self.Tree.SetBranchStatus("Track_z0", 1)
		self.Tree.SetBranchStatus("Track_velX", 1)
		self.Tree.SetBranchStatus("Track_velY", 1)
		self.Tree.SetBranchStatus("Track_velZ", 1)

		self.Tree.GetEntry(self.EventNumber)

		for trackNumber in range(int(self.Tree.NumTracks)):
			x0, y0, z0 = self.Tree.Track_x0[trackNumber], self.Tree.Track_y0[trackNumber], self.Tree.Track_z0[trackNumber]
			vx, vy, vz = self.Tree.Track_velX[trackNumber], self.Tree.Track_velY[trackNumber], self.Tree.Track_velZ[trackNumber]
			xl, yl, zl = det.RecoTrackPoints(x0, y0, z0, vx, vy, vz)
			list_of_trackPt_lists.append( [ physics.RecoTrackPt( xl[n], yl[n], zl[n]  ) for n in range(len(xl)) ] )
			list_of_colors.append(list_of_trackPt_lists[trackNumber][0].c)
		print("digi num hits is ", self.Tree.Digi_numHits)
		for n in range(len(self.Tree.Digi_numHits)):
			x = self.Tree.Digi_x[n]
			y = self.Tree.Digi_y[n]
			z = self.Tree.Digi_z[n]
			self.visEngine.AddPoint( [x, y, z] )
		self.visEngine.TrackDisplay(list_of_trackPt_lists, list_of_colors)
		self.visEngine.Draw()

	def DrawReco_k(self):
		list_of_trackPt_lists = []
		list_of_colors = []
		det = Detector()

		self.Tree.SetBranchStatus("Track_k_x0", 1)
		self.Tree.SetBranchStatus("Track_k_y0", 1)
		self.Tree.SetBranchStatus("Track_k_z0", 1)
		self.Tree.SetBranchStatus("Track_k_velX", 1)
		self.Tree.SetBranchStatus("Track_k_velY", 1)
		self.Tree.SetBranchStatus("Track_k_velZ", 1)

		self.Tree.GetEntry(self.EventNumber)

		for trackNumber in range(len(self.Tree.Track_k_x0)):
			x0, y0, z0 = self.Tree.Track_k_x0[trackNumber], self.Tree.Track_k_y0[trackNumber], self.Tree.Track_k_z0[trackNumber]
			vx, vy, vz = self.Tree.Track_k_velX[trackNumber], self.Tree.Track_k_velY[trackNumber], self.Tree.Track_k_velZ[trackNumber]
			xl, yl, zl = det.RecoTrackPoints(x0, y0, z0, vx, vy, vz)
			list_of_trackPt_lists.append( [ physics.RecoTrackPt( xl[n], yl[n], zl[n]  ) for n in range(len(xl)) ] )
			list_of_colors.append(list_of_trackPt_lists[trackNumber][0].c)
		print("digi num hits is ", self.Tree.Digi_numHits)
		for n in range(len(self.Tree.Digi_numHits)):
			x = self.Tree.Digi_x[n]
			y = self.Tree.Digi_y[n]
			z = self.Tree.Digi_z[n]
			self.visEngine.AddPoint( [x, y, z] )
		self.visEngine.TrackDisplay(list_of_trackPt_lists, list_of_colors)
		self.visEngine.Draw()


	def DrawRecoPoint(self):

		visEngine_local = visualization.Visualizer()

		visEngine_local.vert_vel = self.vert_vel

		visEngine_local.writeDirectory = self.writeDirectory

#		list_of_trackPt_lists = []
#		list_of_colors = []
		det = Detector()

		self.scattered = True

		self.Tree.SetBranchStatus("Track_x0", 1)
		self.Tree.SetBranchStatus("Track_y0", 1)
		self.Tree.SetBranchStatus("Track_z0", 1)
		self.Tree.SetBranchStatus("Track_velX", 1)
		self.Tree.SetBranchStatus("Track_velY", 1)
		self.Tree.SetBranchStatus("Track_velZ", 1)

		self.Tree.SetBranchStatus("Track_k_x0", 1)
		self.Tree.SetBranchStatus("Track_k_y0", 1)
		self.Tree.SetBranchStatus("Track_k_z0", 1)
		self.Tree.SetBranchStatus("Track_k_velX", 1)
		self.Tree.SetBranchStatus("Track_k_velY", 1)
		self.Tree.SetBranchStatus("Track_k_velZ", 1)

		self.Tree.SetBranchStatus("Track_k_hitIndices", 1)
		self.Tree.SetBranchStatus("Track_hitIndices", 1)
		self.Tree.SetBranchStatus("Digi_numHits", 1)
		self.Tree.SetBranchStatus("Digi_x", 1)
		self.Tree.SetBranchStatus("Digi_y", 1)
		self.Tree.SetBranchStatus("Digi_z", 1)

		self.Tree.SetBranchStatus("x_estimates", 1)
		self.Tree.SetBranchStatus("y_estimates", 1)
		self.Tree.SetBranchStatus("z_estimates", 1)

		self.Tree.SetBranchStatus("Hit_particlePx",1)
		self.Tree.SetBranchStatus("Hit_particlePy",1)
		self.Tree.SetBranchStatus("Hit_particlePz",1)

		self.Tree.SetBranchStatus("x_estimates_m", 1)
		self.Tree.SetBranchStatus("y_estimates_m", 1)
		self.Tree.SetBranchStatus("z_estimates_m", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)
		
		if self.vertex:
			self.Tree.SetBranchStatus("Vertex_x",1)
			self.Tree.SetBranchStatus("Vertex_y",1)
			self.Tree.SetBranchStatus("Vertex_z",1)

			self.Tree.SetBranchStatus("Vertex_k_x",1)
			self.Tree.SetBranchStatus("Vertex_k_y",1)
			self.Tree.SetBranchStatus("Vertex_k_z",1)

		self.Tree.GetEntry(self.EventNumber)

		if self.kalman:
			digi_hit_inds = util.unzip(self.Tree.Track_k_hitIndices)

			used_inds = set(self.Tree.Track_k_hitIndices)

		else:
			digi_hit_inds = util.unzip(self.Tree.Track_hitIndices)

			used_inds = set(self.Tree.Track_hitIndices)


#		for trackNumber in range(int(self.Tree.NumTracks)):



		if self.merged:
			numtracks = int(self.Tree.NumTracks_k_m)
			print("number of kalman merged tracks made: ",numtracks)

		elif self.kalman:
			numtracks = len(self.Tree.Track_k_x0)
			print("number of kalman tracks made: ",numtracks)
		else:
			numtracks = int(self.Tree.NumTracks)
			print("number of linear tracks made: ",numtracks)


		if self.merged and self.scattered:
			x_s = util.unzip(self.Tree.x_estimates_m)
			y_s = util.unzip(self.Tree.y_estimates_m)
			z_s = util.unzip(self.Tree.z_estimates_m)

		elif self.kalman and self.scattered:
			x_s = util.unzip(self.Tree.x_estimates)
			y_s = util.unzip(self.Tree.y_estimates)
			z_s = util.unzip(self.Tree.z_estimates)

		for trackNumber in range(numtracks):

#			print(x_s)

			if self.kalman and self.scattered:

				visEngine_local.TrackDisplayPoints(x_s[trackNumber], y_s[trackNumber], z_s[trackNumber])

			else:
				if self.kalman:
					x0, y0, z0 = self.Tree.Track_k_x0[trackNumber], self.Tree.Track_k_y0[trackNumber], self.Tree.Track_k_z0[trackNumber]
					vx, vy, vz = self.Tree.Track_k_velX[trackNumber], self.Tree.Track_k_velY[trackNumber], self.Tree.Track_k_velZ[trackNumber]
				else:
					x0, y0, z0 = self.Tree.Track_x0[trackNumber], self.Tree.Track_y0[trackNumber], self.Tree.Track_z0[trackNumber]
					vx, vy, vz = self.Tree.Track_velX[trackNumber], self.Tree.Track_velY[trackNumber], self.Tree.Track_velZ[trackNumber]

				[xi, yi, zi] = det.FindIntercept(x0, y0, z0, vx, vy, vz) # find intercept with boundary

                #print("the intercepts are x = {}, y = {}, z = {} \n\n".format(xi,yi,zi))

				visEngine_local.TrackDisplayPoints([x0,xi], [y0,yi], [z0,zi])


			if self.used:
        			for ind in digi_hit_inds[trackNumber]:
        				x = self.Tree.Digi_x[ind]
        				y = self.Tree.Digi_y[ind]
        				z = self.Tree.Digi_z[ind]
        				visEngine_local.AddPoint( [x, y, z] )

#			print("hit indices are {}".format(digi_hit_inds))

		num_vertices = 0

		if self.vertex:
			if self.merged:
				num_vertices = len(self.Tree.Vertex_k_m_x)
				description = 'merged kalman'
				
				for n in range(len(self.Tree.Vertex_k_m_x)):
					xv = self.Tree.Vertex_k_m_x[n]
					yv = self.Tree.Vertex_k_m_y[n]
					zv = self.Tree.Vertex_k_m_z[n]
					visEngine_local.AddVertex( [xv, yv, zv] )

			elif self.kalman:
				num_vertices = len(self.Tree.Vertex_k_x)
				description = 'unmerged kalman'

				for n in range(len(self.Tree.Vertex_k_x)):
					xv = self.Tree.Vertex_k_x[n]
					yv = self.Tree.Vertex_k_y[n]
					zv = self.Tree.Vertex_k_z[n]
					visEngine_local.AddVertex( [xv, yv, zv] )

			else:
				num_vertices = len(self.Tree.Vertex_x)
				description = 'linear'

				for n in range(len(self.Tree.Vertex_x)):
					xv = self.Tree.Vertex_x[n]
					yv = self.Tree.Vertex_y[n]
					zv = self.Tree.Vertex_z[n]
					visEngine_local.AddVertex( [xv, yv, zv] )

		if self.unused:
			for n in range(len(self.Tree.Digi_numHits)):
				if not (n in used_inds):
					x = self.Tree.Digi_x[n]
					y = self.Tree.Digi_y[n]
					z = self.Tree.Digi_z[n]
					visEngine_local.AddHit( [x, y, z] )

		print("number of {} vertices made is ".format(description),num_vertices)

		if self.merged:
			visEngine_local.Draw(outname=self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'merged'))
		else:
			visEngine_local.Draw(outname=self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,['linear','kalman'][self.kalman]))



	def RecoVertex(self):

		visEngine_local = visualization.Visualizer()

#		print(self.vert_vel," reco v")
		visEngine_local.vert_vel = self.vert_vel

		visEngine_local.writeDirectory = self.writeDirectory

		self.Tree.SetBranchStatus("Vertex_k_m_x",1)
		self.Tree.SetBranchStatus("Vertex_k_m_y",1)
		self.Tree.SetBranchStatus("Vertex_k_m_z",1)
		self.Tree.SetBranchStatus("Vertex_k_m_trackIndices",1)

		self.Tree.SetBranchStatus("king_move_inds", 1)
		self.Tree.SetBranchStatus("x_estimates_m", 1)
		self.Tree.SetBranchStatus("y_estimates_m", 1)
		self.Tree.SetBranchStatus("z_estimates_m", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)

		self.Tree.SetBranchStatus("vertex_vx_m", 1)
		self.Tree.SetBranchStatus("vertex_vy_m", 1)
		self.Tree.SetBranchStatus("vertex_vz_m", 1)

		self.Tree.SetBranchStatus("Digi_numHits", 1)
		self.Tree.SetBranchStatus("Digi_x", 1)
		self.Tree.SetBranchStatus("Digi_y", 1)
		self.Tree.SetBranchStatus("Digi_z", 1)

		self.Tree.GetEntry(self.EventNumber)

		digi_hit_inds = util.unzip(self.Tree.Track_k_m_hitIndices)
		used_inds = set(self.Tree.Track_k_m_hitIndices)

		used_track_inds = set()

		vert_inds = util.unzip(self.Tree.Vertex_k_m_trackIndices)

#		print(vert_inds)

		colors = ["tab:"+col for col in ['blue','orange','green','red','purple','brown','pink','olive','cyan']]

		x_s = util.unzip(self.Tree.x_estimates_m)
		y_s = util.unzip(self.Tree.y_estimates_m)
		z_s = util.unzip(self.Tree.z_estimates_m)

		vxv = util.unzip(self.Tree.vertex_vx_m) # velocity best estimates at the vertex
		vyv = util.unzip(self.Tree.vertex_vy_m)
		vzv = util.unzip(self.Tree.vertex_vz_m)

#		print("vxv is ", vxv)


		for n in range(len(self.Tree.Vertex_k_m_x)):
			c = n if n < len(colors) else n%len(colors)

			xv = self.Tree.Vertex_k_m_x[n]
			yv = self.Tree.Vertex_k_m_y[n]
			zv = self.Tree.Vertex_k_m_z[n]

			#visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[vxv[n], vyv[n], vzv[n]]} )

			if len(vxv) != 0:
				visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[vxv[n], vyv[n], vzv[n]]} )

			else:
				visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[[]*3]} )

			for trk_ind in vert_inds[n]:
				if self.used:
					for ind in digi_hit_inds[int(trk_ind)]:
						x = self.Tree.Digi_x[ind]
						y = self.Tree.Digi_y[ind]
						z = self.Tree.Digi_z[ind]
						visEngine_local.AddPoint( [[x, y, z], colors[c]] )

				visEngine_local.TrackDisplayPoints(x_s[int(trk_ind)], y_s[int(trk_ind)], z_s[int(trk_ind)], color=colors[c])

				used_track_inds.add(int(trk_ind))


		if self.unused:
			for n in range(len(self.Tree.Digi_numHits)):
				if not (n in used_inds):
					x = self.Tree.Digi_x[n]
					y = self.Tree.Digi_y[n]
					z = self.Tree.Digi_z[n]
					visEngine_local.AddHit( [x, y, z] )

			print(used_track_inds)

			for n in range(self.Tree.NumTracks_k_m):
				if not (n in used_track_inds):
					visEngine_local.TrackDisplayPoints(x_s[n], y_s[n], z_s[n], color='k', opac=0.2)

					print(n)

					for ind in digi_hit_inds[n]:
						x = self.Tree.Digi_x[ind]
						y = self.Tree.Digi_y[ind]
						z = self.Tree.Digi_z[ind]
						visEngine_local.AddPoint( [[x, y, z], 'tab:gray'] )

			king_move_inds = util.unzip(self.Tree.king_move_inds)

			print(king_move_inds)

			for inds in king_move_inds:
				for ind in inds:
					x = self.Tree.Digi_x[ind]
					y = self.Tree.Digi_y[ind]
					z = self.Tree.Digi_z[ind]
					visEngine_local.AddPoint( [[x, y, z], 'm'] )

		visEngine_local.Draw(self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'merged')) #_linear'))


	def VertexAccuracy(self):

		self.Tree.SetBranchStatus("GenParticle_G4index",1)
		self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
		self.Tree.SetBranchStatus("GenParticle_y",1)
		self.Tree.SetBranchStatus("GenParticle_z",1)

		self.Tree.SetBranchStatus("Vertex_k_m_x",1)
		self.Tree.SetBranchStatus("Vertex_k_m_y",1)
		self.Tree.SetBranchStatus("Vertex_k_m_z",1)


		self.Tree.GetEntry(self.EventNumber)

		used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
		dump(used_gens_inds,"used_gens_inds.joblib")
		#used_gens_inds = load("used_gens_inds.joblib")

		#print("Gen indices of used particles ",used_gens_inds)

		vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
			      self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
			      self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10] # only first since all used
									  # particles form one vertex
		#print("vertex truth position is ",vert_truth)

#		print("x zipped is ",self.Tree.Vertex_k_m_x)

#		x = util.unzip(self.Tree.Vertex_k_m_x)
#		y = util.unzip(self.Tree.Vertex_k_m_y)
#		z = util.unzip(self.Tree.Vertex_k_m_z)
		x = self.Tree.Vertex_k_m_x
		y = self.Tree.Vertex_k_m_y
		z = self.Tree.Vertex_k_m_z

#		print("x unzipped is ",x)

		drs = []

		for vert in range(len(x)):
			dx = x[vert] - vert_truth[0]
			dy = y[vert] - vert_truth[1]
			dz = z[vert] - vert_truth[2]

			dr = np.sqrt(dx**2 + dy**2 + dz**2)

			drs.append(dr)

			#print("distance from truth is ",dr)

		return drs

	def PlotMomentum(self, num=100, cut=0):

		#print(self.Tree)

		self.Tree.SetBranchStatus("Hit_particlePx",1)
		self.Tree.SetBranchStatus("Hit_particlePy",1)
		self.Tree.SetBranchStatus("Hit_particlePz",1)
		self.Tree.SetBranchStatus("Hit_particlePdgId", 1)

        # plot momentum distribution
		px,py,pz = [], [], []

		for ev in range(num):
			self.Tree.GetEntry(ev)
			#print(np.array(self.Tree.Hit_particlePx))
			px.extend(np.array(self.Tree.Hit_particlePx))
			py.extend(np.array(self.Tree.Hit_particlePy))
			pz.extend(np.array(self.Tree.Hit_particlePz))
#        		print(px)
#        		print(np.where(px == -1))
#        		for pi in [px, py, pz]:
#        			np.delete(pi, np.where(pi == -1))
		px = np.array(px)
		py = np.array(py)
		pz = np.array(pz)
		p = np.sqrt(px**2 + py**2 + pz**2)

		cap = 1e5

		p[p < cut] = 0
		p[p > cap] = 0

		p = np.delete(p,np.where(p==0))

		#print(p)

		visualization.Histogram(p, Title="Sim Momenta (cut at {} MeV)".format(cut), xaxis="Momentum [MeV]")


	def PlotBeta(self):

		self.Tree.SetBranchStatus("vertex_vx_m", 1)
		self.Tree.SetBranchStatus("vertex_vy_m", 1)
		self.Tree.SetBranchStatus("vertex_vz_m", 1)

		betas = []

		for ev in range(self.Tree.GetEntries()): # event
			self.Tree.GetEntry(ev)

			vx = util.unzip(self.Tree.vertex_vx_m)
			vy = util.unzip(self.Tree.vertex_vy_m)
			vz = util.unzip(self.Tree.vertex_vz_m)

			for v in range(len(vx)): #vertex
				for t in range(len(vx[v])): #track
					betas.append(np.sqrt(vx[v][t]**2 + vy[v][t]**2 + vz[v][t]**2) / physics.c)

		visualization.Histogram(betas, Title="Vertex Velocity Best Estimates", xaxis="beta", log=True, rng=(0,10))


	def PlotVertexBeta(self):

		self.Tree.SetBranchStatus("Vertex_x",1)
		self.Tree.SetBranchStatus("Vertex_y",1)
		self.Tree.SetBranchStatus("Vertex_z",1)
		self.Tree.SetBranchStatus("Vertex_t",1)

		self.Tree.SetBranchStatus("Track_x0", 1)
		self.Tree.SetBranchStatus("Track_y0", 1)
		self.Tree.SetBranchStatus("Track_z0", 1)
		self.Tree.SetBranchStatus("Track_t0", 1)

		self.Tree.SetBranchStatus("Vertex_trackIndices",1)

		betas = []

		for ev in range(self.Tree.GetEntries()): # event
			self.Tree.GetEntry(ev)

			vertex_trackIndices = self.Tree.Vertex_trackIndices
			vertex_trackIndices = util.unzip(vertex_trackIndices)

			#print(vertex_trackIndices)

			if len(vertex_trackIndices) > 0:
				for vert in range(len(vertex_trackIndices)): # vertex
					vert = int(vert)

					vx = self.Tree.Vertex_x[vert]
					vy = self.Tree.Vertex_y[vert]
					vz = self.Tree.Vertex_z[vert]
					vt = self.Tree.Vertex_t[vert]

					for ind in vertex_trackIndices[vert]: # track
						x = self.Tree.Track_x0[ind]
						y = self.Tree.Track_y0[ind]
						z = self.Tree.Track_z0[ind]
						t = self.Tree.Track_t0[ind]

						beta = np.sqrt(((x-vx)**2 + (y-vy)**2 + (z-vz)**2) / (physics.c**2 * (t-vt)**2))

						betas.append(beta)

		visualization.Histogram(betas)


	def PlotHitAccuracy(self, num=100):

		self.Tree.SetBranchStatus("track_pdgs", 1)
		self.Tree.SetBranchStatus("track_ids", 1)

		dist_of_particles = []
		dist_of_ids = []

		for ev in range(num):
			self.Tree.GetEntry(ev)

			ids_list = util.unzip(self.Tree.track_ids)
			pdgs_list = util.unzip(self.Tree.track_pdgs)

			for track_pdg in pdgs_list:
				dist_of_particles.append(len(set(track_pdg)))

			for track_ids in ids_list:
				dist_of_ids.append(len(set(track_ids)))

		visualization.Histogram(dist_of_particles,Title="Number of PDG IDs per track",xaxis="Number of IDs")
		visualization.Histogram(dist_of_ids,Title="Number of Track IDs per track",xaxis="Number of IDs")



#			numtracks = len(self.Tree.Track_k_x0)
#
#			for trackNumber in range(numtracks):
#				pdgs = set()
#
#				for hit in


	def find_with_bool(self, bool_func, Op=0):

		inds = []

		for ev in range(int(self.Tree.GetEntries())):
			self.Tree.GetEntry(ev)

			if bool_func(self.Tree, op=Op):
				inds.append(ev)

		return inds


	def find_missed_tracks(self, num, mult):

		self.Tree.SetBranchStatus("NumTracks", 1)
		self.Tree.SetBranchStatus("NumTracks_k", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)
		missed_inds_k = []
		missed_inds_lin = []
		high_mult = []

		for ev in range(num[0],num[1]):

			self.Tree.GetEntry(ev)

			if self.merged:
				numtracks_k = self.Tree.NumTracks_k_m
			else:
				numtracks_k = self.Tree.NumTracks_k

			numtracks = self.Tree.NumTracks

			if numtracks != 0 and numtracks_k == 0:

				missed_inds_k.append(ev)

			if numtracks == 0 and numtracks_k != 0:

				missed_inds_lin.append(ev)

			if numtracks_k > mult:

				high_mult.append(ev)


#		print("kalman events with multiplicity > {}\n".format(mult),high_mult)

		return missed_inds_k, missed_inds_lin, high_mult

	def find_bracket(self, num):

		self.Tree.SetBranchStatus("NumTracks", 1)
		self.Tree.SetBranchStatus("NumTracks_k", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)

		bracket = []

		for ev in range(num[0],num[1]):

			self.Tree.GetEntry(ev)

			numtracks_k_m = int(self.Tree.NumTracks_k_m)
			numtracks_k = int(self.Tree.NumTracks_k)

			# if numtracks_k >= 2 and numtracks_k_m < 2:
			# 	bracket.append(ev)

			if numtracks_k_m > 1:
				bracket.append(ev)

		return bracket


	def GetRecoInfo(self):
		self.Tree.SetBranchStatus("Digi_x", 1)
		self.Tree.SetBranchStatus("Digi_y", 1)
		self.Tree.SetBranchStatus("Digi_z", 1)
		self.Tree.SetBranchStatus("NumTracks", 1)
		self.Tree.SetBranchStatus("NumVertices", 1)
		self.Tree.SetBranchStatus("Vertex_x", 1)
		self.Tree.SetBranchStatus("Vertex_y", 1)
		self.Tree.SetBranchStatus("Vertex_z", 1)
		self.Tree.SetBranchStatus("Vertex_t", 1)
		self.Tree.SetBranchStatus("Vertex_ErrorY", 1)
		self.Tree.SetBranchStatus("Vertex_ErrorT", 1)
		self.Tree.SetBranchStatus("Track_velX", 1)
		self.Tree.SetBranchStatus("Track_velY", 1)
		self.Tree.SetBranchStatus("Track_velZ", 1)
		self.Tree.SetBranchStatus("Track_x0", 1)
		self.Tree.SetBranchStatus("Track_y0", 1)
		self.Tree.SetBranchStatus("Track_z0", 1)
		self.Tree.SetBranchStatus("Track_t0", 1)
		self.Tree.SetBranchStatus("Track_missingHitLayer", 1)
		self.Tree.SetBranchStatus("Track_expectedHitLayer", 1)
		self.Tree.SetBranchStatus("track_ipDistance", 1)
		self.Tree.SetBranchStatus("Track_hitIndices", 1)
		self.Tree.SetBranchStatus("Track_beta", 1)
		self.Tree.SetBranchStatus("Track_ErrorBeta", 1)

		self.Tree.GetEntry(self.EventNumber)

		print("Number of Tracks: " + str(self.Tree.NumTracks))

		associated_digis = util.unzip(self.Tree.Track_hitIndices)
		missing_hits     = util.unzip(self.Tree.Track_missingHitLayer)
		expected_hits    = util.unzip(self.Tree.Track_expectedHitLayer)

		for n in range(int(self.Tree.NumTracks)):
			print("**Track: " + str(n) + "**")
			print("Start Point: (" + str(self.Tree.Track_x0[n]) + ", " + str(self.Tree.Track_y0[n]) + ", " + str(self.Tree.Track_z0[n]) + ")")
			print("Velocity:    (" + str(self.Tree.Track_velX[n]) + ", " + str(self.Tree.Track_velY[n]) + ", " + str(self.Tree.Track_velZ[n]) + ")")
			print("Beta: " + str(self.Tree.Track_beta[n]) + " +/- " + str(self.Tree.Track_ErrorBeta[n]))
			print("Digis: ")
			for digi_index in associated_digis[n]:
				print("--Digi " + str(digi_index))
				print("--(" + str(self.Tree.Digi_x[digi_index]) + ", " + str(self.Tree.Digi_y[digi_index]) + ", " + str(self.Tree.Digi_z[digi_index]) + ")" )
#			print("Missing Hits in Layers:  " + str(missing_hits[n]))
#			print("Expected Hits in Layers: " + str(expected_hits[n]))



	def TrackEfficiency(self):
		''' plot distribution of reconstructed lowest hit location'''

		self.Tree.SetBranchStatus("Track_k_m_x0", 1)
		self.Tree.SetBranchStatus("Track_k_m_z0", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)

		x0, z0 = [], []

		has_a_track = 0

		for ev in range(self.Tree.GetEntries()): # event
			self.Tree.GetEntry(ev)

			for tr in range(self.Tree.NumTracks_k_m):
				x0.append(self.Tree.Track_k_m_x0[tr])
				z0.append(self.Tree.Track_k_m_z0[tr])

			if self.Tree.NumTracks_k_m > 0:
				has_a_track += 1

		#plt.hist2d(trackx0, trackz0)

		hist, bin_x, bin_y = np.histogram2d(x0,z0,bins=150)

		fig, ax = plt.subplots(figsize=(8,5))
		plt.imshow(hist,alpha=0.5)
		plt.colorbar(orientation='vertical')
		plt.savefig("Efficiency.png")

		efficiency = has_a_track / self.Tree.GetEntries()

		print("efficiency is ", efficiency)


	def TrackResolution(self):
		''' plot distribution of truth vs reconstructed for angle and lowest hit location '''

		self.Tree.SetBranchStatus("Track_k_m_x0", 1)
		self.Tree.SetBranchStatus("Track_k_m_z0", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)



