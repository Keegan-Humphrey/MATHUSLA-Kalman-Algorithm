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

		self.Tree.SetBranchStatus("GenParticle_G4index",1)
		self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
		self.Tree.SetBranchStatus("GenParticle_y",1)
		self.Tree.SetBranchStatus("GenParticle_z",1)
		self.Tree.SetBranchStatus("Vertex_x", 1)
		self.Tree.SetBranchStatus("Vertex_y", 1)
		self.Tree.SetBranchStatus("Vertex_z", 1)

		visEngine_local = visualization.Visualizer()

		visEngine_local.writeDirectory = self.writeDirectory

		used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
		dump(used_gens_inds,"used_gens_inds.joblib")
                #used_gens_inds = load("used_gens_inds.joblib")

		vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
                              self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
                              self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10] # only first since all used

                #truth vertex location
		visEngine_local.AddVertex( {'point':[vert_truth[0], vert_truth[1], vert_truth[2]], 'col':'tab:green', 'vert vel':[[],[],[]]} )

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


		'''
		for n in range(len(self.Tree.Vertex_x)):
			xv = self.Tree.Vertex_x[n]
			yv = self.Tree.Vertex_y[n]
			zv = self.Tree.Vertex_z[n]

			# made vertex passed to visualiser
			visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':'tab:red', 'vert vel':[[]*3]} )
		'''

		self.num_truth = num_truth
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


	def find_with_bool(self, bool_func, Op=0):
		''' find indices of events that satisfy a criteria given by bool_func (boolean valued function)'''

		inds = []

		for ev in range(int(self.Tree.GetEntries())):
			self.Tree.GetEntry(ev)

			if bool_func(self.Tree, op=Op):
				inds.append(ev)

		return inds


	def RecoVertex(self):
		''' visualise tracks and vertices made by kalman filter, and the digi hits and truth decay locations '''

		visEngine_local = visualization.Visualizer()

		visEngine_local.vert_vel = self.vert_vel

		visEngine_local.writeDirectory = self.writeDirectory

		self.Tree.SetBranchStatus("GenParticle_G4index",1)
		self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
		self.Tree.SetBranchStatus("GenParticle_y",1)
		self.Tree.SetBranchStatus("GenParticle_z",1)

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

		# indices of used digi hits (organised by track and organised into a set)
		digi_hit_inds = util.unzip(self.Tree.Track_k_m_hitIndices)
		used_inds = set(self.Tree.Track_k_m_hitIndices)

		used_track_inds = set()

		# indices of tracks used in the made vertices
		vert_inds = util.unzip(self.Tree.Vertex_k_m_trackIndices)

#		print(vert_inds)

		# track best estimates
		x_s = util.unzip(self.Tree.x_estimates_m)
		y_s = util.unzip(self.Tree.y_estimates_m)
		z_s = util.unzip(self.Tree.z_estimates_m)

		# velocity best estimates at the vertex
		vxv = util.unzip(self.Tree.vertex_vx_m)
		vyv = util.unzip(self.Tree.vertex_vy_m)
		vzv = util.unzip(self.Tree.vertex_vz_m)

		# get truth vertex info (decay location of parent particle used in the event)
		used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
		#dump(used_gens_inds,"used_gens_inds.joblib")
		#used_gens_inds = load("used_gens_inds.joblib")

		# Vertex truth location (translated from Pythia coordinate system)
		vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
			      self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
			      self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10]

		# pass truth vertex to visualiser
		visEngine_local.AddVertex( {'point':[vert_truth[0], vert_truth[1], vert_truth[2]], 'col':'tab:green', 'vert vel':[[],[],[]]} )

		# track and vertex colors
		colors = ["tab:"+col for col in ['blue','orange','red','purple','brown','pink','olive','cyan']]

		# pass vertices, tracks, and hits used in them to visualiser
		# (use one color per vertex and its tracks and their hits)
		for n in range(len(self.Tree.Vertex_k_m_x)):
			# choose index of a color to use for visualiser
			c = n if n < len(colors) else n%len(colors)

			xv = self.Tree.Vertex_k_m_x[n]
			yv = self.Tree.Vertex_k_m_y[n]
			zv = self.Tree.Vertex_k_m_z[n]

			# made vertex passed to visualiser
			if len(vxv) != 0:
				visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[vxv[n], vyv[n], vzv[n]]} )

			else:
				visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[[]*3]} )

			# loop over tracks in the current vertex
			for trk_ind in vert_inds[n]:
				if self.used:
					for ind in digi_hit_inds[int(trk_ind)]:
						x = self.Tree.Digi_x[ind]
						y = self.Tree.Digi_y[ind]
						z = self.Tree.Digi_z[ind]

						# pass each hit used in the current track to visualiser
						visEngine_local.AddPoint( [[x, y, z], colors[c]] )

				# pass track to visualiser
				visEngine_local.TrackDisplayPoints(x_s[int(trk_ind)], y_s[int(trk_ind)], z_s[int(trk_ind)], color=colors[c])

				used_track_inds.add(int(trk_ind))


		if self.unused:
			# pass hits that weren't used in a track to visualiser
			for n in range(len(self.Tree.Digi_numHits)):
				if not (n in used_inds):
					x = self.Tree.Digi_x[n]
					y = self.Tree.Digi_y[n]
					z = self.Tree.Digi_z[n]
					visEngine_local.AddHit( [x, y, z] )

			#print(used_track_inds)

			# pass tracks not used in a vertex and their used hits to visualiser
			for n in range(self.Tree.NumTracks_k_m):
				if not (n in used_track_inds):
					visEngine_local.TrackDisplayPoints(x_s[n], y_s[n], z_s[n], color='k', opac=0.2)

					#print(n)

					for ind in digi_hit_inds[n]:
						x = self.Tree.Digi_x[ind]
						y = self.Tree.Digi_y[ind]
						z = self.Tree.Digi_z[ind]
						visEngine_local.AddPoint( [[x, y, z], 'tab:gray'] )

			king_move_inds = util.unzip(self.Tree.king_move_inds)

			print('king move indices are ',king_move_inds)

			# show hits dropped from the hit pool (via king moves algorithm)
			for inds in king_move_inds:
				for ind in inds:
					x = self.Tree.Digi_x[ind]
					y = self.Tree.Digi_y[ind]
					z = self.Tree.Digi_z[ind]
					visEngine_local.AddPoint( [[x, y, z], 'm'] )

		# draw the plot and save (show) it
		visEngine_local.Draw(self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'merged')) #_linear'))



	def VertexAccuracy(self):
		''' get differences in each coordinate between decay location (truth vertex) and the closest made vertex'''

		self.Tree.SetBranchStatus("GenParticle_G4index",1)
		self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
		self.Tree.SetBranchStatus("GenParticle_y",1)
		self.Tree.SetBranchStatus("GenParticle_z",1)

		self.Tree.SetBranchStatus("Vertex_k_m_x",1)
		self.Tree.SetBranchStatus("Vertex_k_m_y",1)
		self.Tree.SetBranchStatus("Vertex_k_m_z",1)
		self.Tree.SetBranchStatus("Vertex_k_m_ErrorX",1)
		self.Tree.SetBranchStatus("Vertex_k_m_ErrorY",1)
		self.Tree.SetBranchStatus("Vertex_k_m_ErrorZ",1)

		dx, dy, dz = [], [], []
		ex, ey, ez = [], [], []
		px, py, pz = [], [], []
		chi = []

		max = 1e3

		res_info_txt = []

		for ev_num in range(self.Tree.GetEntries()):

			print("event {}".format(ev_num)) if ev_num % 100 == 0 else None

			self.Tree.GetEntry(ev_num)

			used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
			#dump(used_gens_inds,"used_gens_inds.joblib")
			#used_gens_inds = load("used_gens_inds.joblib")

			# particles form one vertex only first since all used
			vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
				      self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
				      self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10]

			x = self.Tree.Vertex_k_m_x
			y = self.Tree.Vertex_k_m_y
			z = self.Tree.Vertex_k_m_z
			_ex = self.Tree.Vertex_k_m_ErrorX
			_ey = self.Tree.Vertex_k_m_ErrorY
			_ez = self.Tree.Vertex_k_m_ErrorZ

			min_dist = 1e6
			min_index = -1

			for vert in range(len(x)):
				_dx = x[vert] - vert_truth[0]
				_dy = y[vert] - vert_truth[1]
				_dz = z[vert] - vert_truth[2]

				dr = np.sqrt(_dx**2 + _dy**2 + _dz**2)

				if dr < min_dist:
					min_index = vert

			if min_index != -1:
				dx.append(x[min_index] - vert_truth[0])
				dy.append(y[min_index] - vert_truth[1])
				dz.append(z[min_index] - vert_truth[2])
				ex.append(_ex[min_index])
				ey.append(_ey[min_index])
				ez.append(_ez[min_index])
				px.append(dx[-1] / ex[-1])
				py.append(dy[-1] / ey[-1])
				pz.append(dz[-1] / ez[-1])
				chi.append(dx[-1]**2 / ex[-1]**2 + dy[-1]**2 / ey[-1]**2 + dz[-1]**2 / ez[-1]**2)

				res_info_txt.append("\nev num is " + str(ev_num)+ '\n')
				res_info_txt.append("chi is " + str(chi[-1]) + '\n')
				res_info_txt.append("pull X is " + str(px[-1]) + '\n')
				res_info_txt.append("pull Y is " + str(py[-1]) + '\n')
				res_info_txt.append("pull Z is " + str(pz[-1]) + '\n')
				res_info_txt.append("error X is " + str(ex[-1]) + '\n')
				res_info_txt.append("error Y is " + str(ey[-1]) + '\n')
				res_info_txt.append("error Z is " + str(ez[-1]) + '\n')
				res_info_txt.append("diff X is " + str(dx[-1]) + '\n')
				res_info_txt.append("diff Y is " + str(dy[-1]) + '\n')
				res_info_txt.append("diff Z is " + str(dz[-1]) + '\n')

		f = open("res_info.txt","w+")
		f.writelines(res_info_txt)
		f.close()

		vis_engine = visualization.Histogram(dx, rng=(-max,max), Title='Vertex X Resolution', \
                	xaxis='Vertex X - Decay Location X [cm]', fname='resolutionX.png')
		vis_engine = visualization.Histogram(dy, rng=(-max,max), Title='Vertex Y Resolution', \
                	xaxis='Vertex Y - Decay Location Y [cm]', fname='resolutionY.png')
		vis_engine = visualization.Histogram(dz, rng=(-max,max), Title='Vertex Z Resolution', \
                	xaxis='Vertex Z - Decay Location Z [cm]', fname='resolutionZ.png')
		vis_engine = visualization.Histogram(ex, rng=(0,max), Title='Vertex X Error', \
                	xaxis='Error X [cm]', fname='errorX.png')
		vis_engine = visualization.Histogram(ey, rng=(0,max), Title='Vertex Y Error', \
                	xaxis='Error Y [cm]', fname='errorY.png')
		vis_engine = visualization.Histogram(ez, rng=(0,max), Title='Vertex Z Error', \
                	xaxis='Error Z [cm]', fname='errorZ.png')
		vis_engine = visualization.Histogram(px, rng=(-100,100), Title='Vertex X Pull', \
                	xaxis='Pull X []', fname='pullX.png')
		vis_engine = visualization.Histogram(py, rng=(-100,100), Title='Vertex Y Pull', \
                	xaxis='Pull Y []', fname='pullY.png')
		vis_engine = visualization.Histogram(pz, rng=(-100,100), Title='Vertex Z Pull', \
                	xaxis='Pull Z []', fname='pullZ.png')
		vis_engine = visualization.Histogram(chi, rng=(0,200), Title='Vertex chi', \
                	xaxis='Chi []', fname='chi.png')



	def PlotMomentum(self, num=100, cut=0):
		''' plot momentum distribution of sim hits '''

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

		visualization.Histogram(p, Title="Sim Momenta (cut at {} MeV)".format(cut), xaxis="Momentum [MeV]")


	def PlotBeta(self):
		''' plot vertex velocity best estimates (only if kalman vertexer is used)'''

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
		''' Plot beta required for a particle to move from bottom of a track to location of the vertex it belongs to'''

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
		self.Tree.SetBranchStatus("Track_k_m_y0", 1)
		self.Tree.SetBranchStatus("Track_k_m_z0", 1)

		self.Tree.SetBranchStatus("Track_k_velX", 1)
		self.Tree.SetBranchStatus("Track_k_velY", 1)
		self.Tree.SetBranchStatus("Track_k_velZ", 1)
		self.Tree.SetBranchStatus("NumTracks_k_m", 1)

		self.Tree.SetBranchStatus("Hit_x", 1)
		self.Tree.SetBranchStatus("Hit_y", 1)
		self.Tree.SetBranchStatus("Hit_z", 1)

		self.Tree.SetBranchStatus("Hit_particlePx", 1)
		self.Tree.SetBranchStatus("Hit_particlePy", 1)
		self.Tree.SetBranchStatus("Hit_particlePz", 1)

		min_distances = []
		angle_difference = []

		for ev in range(self.Tree.GetEntries()): # event
			self.Tree.GetEntry(ev)

			for tr in range(self.Tree.NumTracks_k_m):
				min_sim_distance = 1e6
				min_sim_index = -1

				for sim in range(len(self.Tree.Hit_x)):
					distance = (self.Tree.Track_k_m_x0[tr] - self.Tree.Hit_x[sim])**2 \
						 + (self.Tree.Track_k_m_y0[tr] - self.Tree.Hit_y[sim])**2 \
						 + (self.Tree.Track_k_m_z0[tr] - self.Tree.Hit_z[sim])**2

					if np.sqrt(distance) < min_sim_distance:
						min_sim_distance = np.sqrt(distance)
						min_sim_index = sim

				min_distances.append(min_sim_distance)

				momentum = np.array([self.Tree.Hit_particlePx[min_sim_index],
						     self.Tree.Hit_particlePy[min_sim_index],
						     self.Tree.Hit_particlePz[min_sim_index]])

				track_vel = np.array([self.Tree.Track_k_velX[tr],
						      self.Tree.Track_k_velY[tr],
						      self.Tree.Track_k_velZ[tr]])

				angle = np.sum(momentum * track_vel) \
					/ np.sqrt(np.sum(track_vel**2) * np.sum(momentum**2))

				angle_difference.append(angle)

		visualization.Histogram(min_distances,Title='Truth vs Track distance',xaxis='Distance [cm]',fname='distance.png')
		visualization.Histogram(angle_difference,Title='Truth vs Track angle',xaxis='cos(angle) []',log=True,fname='angle.png')




