from detector import Detector 
import os
import ROOT as root
import util
from event import Event
import numpy as np

class H_mumu_Analyzer:

	plot_dir = ""
	NCUTS = 6
	events_passing_cuts = [0.0 for n in range(NCUTS+1)]
	det = Detector()

	def __init__(self, loop_dir):
		#self.files = [loop_dir] # global path and filename
		self.files = util.GetFilesInDir(loop_dir)
		self.passed_files = []
		self.passed_events = []
		self.floor_hit_location = root.TH2D("floor_hit_location", "floor hit x,z", 1000, -5100., 5100., 1000, 6900., 17100. )

	def SetPlotDir(self, dirname):
		self.plot_dir = dirname

	def Analyze(self, tree_name="integral_tree"):
		self.tree_name = tree_name
		for file in self.files:
			passed_events = []
			print("Working in file: " + file)
			tfile = root.TFile.Open(file)
			self.events_passing_cuts_byfile = [0.0 for n in range(self.NCUTS+1)]
			print(tree_name)
			print(tfile.Get(tree_name))
			print("tfile is ",tfile)
			self.InitTree(tfile.Get(tree_name))
			for eventNumber in range(self.Tree.GetEntries()):
				self.Tree.GetEntry(eventNumber)
				if self.Selection():
					print(file + ": " + str(eventNumber))
					passed_events.append(eventNumber)
			if len(passed_events) > 0:
				self.passed_files.append(file)
				self.passed_events.append(passed_events)

			print(self.events_passing_cuts_byfile)

		print(self.passed_events)
		print("H_mumu Analyzer Results:")
		print(self.events_passing_cuts)
		c1 = root.TCanvas("c1")
		self.floor_hit_location.SetMarkerSize(2)
		self.floor_hit_location.SetMarkerStyle(6)
		self.floor_hit_location.Draw()
		c1.Print("floor_hit.png")

	def StudyPassedEvents(self, n):
		file = self.passed_files[n]
		tfile = root.TFile.Open(file)
		self.Tree = tfile.Get(self.tree_name)
		for eventNumber in self.passed_events[n]:
			self.Tree.GetEntry(eventNumber)
			currEvent = Event(self.Tree, eventNumber)
			#currEvent.ExtractTruthPhysics()
			#currEvent.Print()
			currEvent.GetRecoInfo()
			currEvent.DrawReco()

	def InitTree(self, tree):
		self.Tree = tree
		self.Tree.SetBranchStatus("*", 0)
		self.Tree.SetBranchStatus("Digi_x", 1)
		self.Tree.SetBranchStatus("Digi_y", 1)
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
		self.Tree.SetBranchStatus("Track_ErrorY0", 1)
		self.Tree.SetBranchStatus("Track_ErrorT0", 1)
		self.Tree.SetBranchStatus("Track_missingHitLayer", 1)
		self.Tree.SetBranchStatus("Track_expectedHitLayer", 1)
		self.Tree.SetBranchStatus("track_ipDistance", 1)
		self.Tree.SetBranchStatus("Track_hitIndices", 1)
		self.Tree.SetBranchStatus("GenParticle_G4index",1)
		self.Tree.SetBranchStatus("GenParticle_pdgid",1)
		self.Tree.SetBranchStatus("GenParticle_x",1)
		self.Tree.SetBranchStatus("GenParticle_y",1)
		self.Tree.SetBranchStatus("GenParticle_z",1)
		self.Tree.SetBranchStatus("GenParticle_time",1)


	def Trigger(self):
		if len(self.Tree.Digi_x) < 3:
			return False 
		return True

	def Selection(self):
		
		if not self.Trigger():
			return False

		###########################################
		#counting total events
		self.events_passing_cuts[0] += 1.0
		self.events_passing_cuts_byfile[0] += 1.0
		###########################################
		
		###########################################
		#ntracks cut
		if (self.Tree.NumTracks < 2):
			return False

		self.events_passing_cuts[1] += 1.0
		self.events_passing_cuts_byfile[1] += 1.0
		###########################################

		###########################################
		#floor veto w/ expected hit cuts
		for hity in self.Tree.Digi_y:
			if self.det.inLayer(hity) < 2:
				return False

		expected_hits = util.unzip(self.Tree.Track_expectedHitLayer)

		bottom_layer_expected_hits = []

		for exp_list in expected_hits:
			for val in exp_list:
				if val < 2:
					bottom_layer_expected_hits.append(val)

		if len(bottom_layer_expected_hits) < 3:
			return False


		self.events_passing_cuts[2] += 1.0
		self.events_passing_cuts_byfile[2] += 1.0
		###########################################


		####
		####

		x00, y00, z00 = self.Tree.Track_x0[0], self.Tree.Track_y0[0], self.Tree.Track_z0[0]
		x01, y01, z01 = self.Tree.Track_x0[1], self.Tree.Track_y0[1], self.Tree.Track_z0[1]

		vx0, vy0, vz0 = self.Tree.Track_velX[0], self.Tree.Track_velY[0], self.Tree.Track_velZ[0]
		vx1, vy1, vz1 = self.Tree.Track_velX[1], self.Tree.Track_velY[1], self.Tree.Track_velZ[1]

		floor_y = 6002.5

		delt0 = (y00 - floor_y)/vy0
		delt1 = (y01 - floor_y)/vy1

		expected_x0 = x00 + delt0*vx0
		expected_x1 = x01 + delt1*vx1
		expected_z0 = z00 + delt0*vz0
		expected_z1 = z01 + delt1*vz1

		#plotting the location of these hits
		self.floor_hit_location.Fill(expected_x0, expected_z0)
		self.floor_hit_location.Fill(expected_x1, expected_z1)


		####
		####

		###########################################
		#nvertices cut
		if self.Tree.NumVertices == 0:
			return False

		self.events_passing_cuts[3] += 1.0
		self.events_passing_cuts_byfile[3] += 1.0
		###########################################

		###########################################
		#fiducial vertex cut
		if not self.det.inBox(self.Tree.Vertex_x[0], self.Tree.Vertex_y[0], self.Tree.Vertex_z[0]):
			return False

		self.events_passing_cuts[4] += 1.0
		self.events_passing_cuts_byfile[4] += 1.0
		###########################################

		

		###########################################
		#vertex before track cut

		vtxTrackConsistencyY = max( [ (self.Tree.Vertex_y[0] - self.Tree.Track_y0[n])/self.Tree.Track_ErrorY0[n] for n in range(int(self.Tree.NumTracks)) ] )
		#vtxTrackConsistencyT = max( [ (self.Tree.Vertex_t[0] - self.Tree.Track_t0[n])/self.Tree.Track_ErrorT0[n] for n in range(int(self.Tree.NumTracks)) ] )

		if vtxTrackConsistencyY > 1.0:
			return

		self.events_passing_cuts[5] += 1.0
		self.events_passing_cuts_byfile[5] += 1.0
		###########################################

		###########################################
		#missing hits in upper layers

		trackn = 0
		vertex_first_layer = self.det.nextLayer(self.Tree.Vertex_y[0])
		for layern in self.Tree.Track_missingHitLayer:
			if layern >= vertex_first_layer:
				return False

		self.events_passing_cuts[6] += 1.0
		self.events_passing_cuts_byfile[6] += 1.0

		#note the cut below isnt necessary when requiring no missing hits
		###########################################

		###########################################
		#tracks in vertex start in same layer

		#track_hit_yvals = [ [] for i in range(len(self.Tree.Track_x0))]
		#trackn = 0
		#for hitn in self.Tree.Track_hitIndices:
		#	if hitn == -1:
		#		trackn += 1
		#	else:
		#		track_hit_yvals[trackn].append(self.Tree.Digi_y[hitn])

		#min_layers = [ self.det.inLayer(min(yvals_list)) for yvals_list in track_hit_yvals ]

		#veto = False

		#start = min_layers[0]

		#for minval in min_layers:
		#	if not minval==start:
		#		#check if there is expected hit in that layer
		#		return False

		#self.events_passing_cuts[7] += 1.0
		#self.events_passing_cuts_byfile[] += 1.0
		###########################################


		return True


	def Plot1D(self, branch_name):
		plotvar = 0.0
		self.Tree.SetBranchStatus("")

	def PlotSelection(self):
		if (int(self.Tree.NumVertices) == 1 and int(self.Tree.NumTracks) == 2):
			if self.det.inBox(self.Tree.Vertex_x[0], self.Tree.Vertex_y[0], self.Tree.Vertex_z[0]):
				return True
		return False

	def Plot(self, tree_name="integral_tree"):
		plotx = root.TH1D("x_res", "Vertex X Resolution (truth-actual)", 100, -500, 500)
		ploty = root.TH1D("y_res", "Vertex Y Resolution (truth-actual)", 100, -500, 500)
		plotz = root.TH1D("z_res", "Vertex Z Resolution (truth-actual)", 100, -500, 500)
		plott = root.TH1D("t_res", "Vertex t Resolution (truth-actual)", 100, -500, 500)

		self.tree_name = tree_name
		for file in self.files:
			print("Working in file: " + file)
			tfile = root.TFile.Open(file)
			self.InitTree(tfile.Get(tree_name))
			for eventNumber in range(self.Tree.GetEntries()):
				self.Tree.GetEntry(eventNumber)
				if self.PlotSelection():
					vx = self.Tree.Vertex_x[0]
					#evx = self.Tree.Vertex_ErrorX[0]
					vy = self.Tree.Vertex_y[0]
					#evy = self.Tree.Vertex_ErrorY[0]
					vz = self.Tree.Vertex_z[0]
					#evz = self.Tree.Vertex_ErrorZ[0]
					vt = self.Tree.Vertex_t[0]
					#evt = self.Tree.Vertex_ErrorT[0]
					for k in range(int(len(self.Tree.GenParticle_G4index))):
						if (self.Tree.GenParticle_G4index[k] == 1) and int(np.absolute(self.Tree.GenParticle_pdgid[k] == 13)) :
							gen_x = self.Tree.GenParticle_y[k]/10.
							plotx.Fill( (gen_x-vx))#/evx )
							gen_y = self.Tree.GenParticle_x[k]/10.
							ploty.Fill( (gen_y-vy))#/evy )
							gen_z = self.Tree.GenParticle_z[k]/10.
							plotz.Fill( (gen_z-vz))#/evz )
							gen_t = self.Tree.GenParticle_time[k]
							plott.Fill( (gen_t-vt))#/evt )

							print( [(vt-gen_t), (vx-gen_x), (vy-gen_y), (vz-gen_z)])
							break
					

		c1 = root.TCanvas("c1")
		plotx.Draw()
		plotx.GetXaxis().SetTitle("distance [cm]")
		c1.Print("plotxcm.png")

		c2 = root.TCanvas("c2")
		ploty.Draw()
		ploty.GetXaxis().SetTitle("distance [cm]")
		c2.Print("plotycm.png")

		c3 = root.TCanvas("c3")
		plotz.Draw()
		plotz.GetXaxis().SetTitle("distance [cm]")
		c3.Print("plotzcm.png")

		c4 = root.TCanvas("c4")
		plott.Draw()
		plott.GetXaxis().SetTitle("time difference [ns]")
		c4.Print("plottcm.png")


					

###############################################################################################################################################################################
###############################################################################################################################################################################


class K_Long_Analayzer:

	plot_dir = ""
	NCUTS = 4
	events_passing_cuts = [0.0 for n in range(NCUTS+1)]
	det = Detector()

	def __init__(self, loop_dir):
		self.files = util.GetFilesInDir(loop_dir)
		self.passed_files = []
		self.passed_events= []

	def SetPlotDir(self, dirname):
		self.plot_dir = dirname

	def Plot(self, tree_name="integral_tree"):
		self.tree_name = tree_name
		plot = root.TH1D("plot", "Track Beta for events w/ vertex", 20, 0.60, 1.20)
		for file in self.files:
			print("Working in file: " + file)
			tfile = root.TFile.Open(file)
			self.InitTree(tfile.Get(tree_name))
			for eventNumber in range(self.Tree.GetEntries()):
				self.Tree.GetEntry(eventNumber)
				plotif, val = self.SelectionForPlot()
				if plotif:
					plot.Fill(val)
		c1 = root.TCanvas("c1")
		plot.Draw()
		c1.Print(self.plot_dir + "plot1.png", ".png")

	def SelectionForPlot(self):
		if not self.Trigger():
			return False, 0

		###########################################
		#counting total events
		
		###########################################
		#ntracks cut
		if (self.Tree.NumTracks < 2):
			return False, 0


		###########################################
		#nvertices cut
		if self.Tree.NumVertices == 0:
			return False, 0


		###########################################
		#fiducial vertex cut
		if not self.det.inBox(self.Tree.Vertex_x[0], self.Tree.Vertex_y[0], self.Tree.Vertex_z[0]):
			return False, 0


		###########################################
		#floor veto w/ expected hit cuts
		for hity in self.Tree.Digi_y:
			if self.det.inLayer(hity) < 2:
				return False, 0

		expected_hits = util.unzip(self.Tree.Track_expectedHitLayer)

		bottom_layer_expected_hits = []

		for exp_list in expected_hits:
			for val in exp_list:
				if val < 2:
					bottom_layer_expected_hits.append(val)

		if len(bottom_layer_expected_hits) < 3:
			return False, 0


		###########################################
		#vertex before track cut

		
		return True, min(self.Tree.Track_beta)


	def Analyze(self, tree_name="integral_tree"):
		self.tree_name = tree_name
		for file in self.files:
			passed_events = []
			print("Working in file: " + file)
			tfile = root.TFile.Open(file)
			self.InitTree(tfile.Get(tree_name))
			for eventNumber in range(self.Tree.GetEntries()):
				self.Tree.GetEntry(eventNumber)
				if self.Selection():
					print(file + ": " + str(eventNumber))
					passed_events.append(eventNumber)
			if len(passed_events) > 0:
				self.passed_files.append(file)
				self.passed_events.append(passed_events)

			print(self.events_passing_cuts_byfile)
		
		print(self.passed_events)
		print("H_mumu Analyzer Results:")
		print(self.events_passing_cuts)

	def StudyPassedEvents(self, n):
		file = self.passed_files[n]
		tfile = root.TFile.Open(file)
		self.Tree = tfile.Get(self.tree_name)
		for eventNumber in self.passed_events[n]:
			self.Tree.GetEntry(eventNumber)
			currEvent = Event(self.Tree, eventNumber)
			#currEvent.ExtractTruthPhysics()
			#currEvent.Print()
			currEvent.GetRecoInfo()
			currEvent.DrawReco()
		

	def InitTree(self, tree):
		self.Tree = tree
		#self.Tree.SetBranchStatus("*", 0)
		self.Tree.SetBranchStatus("Digi_x", 1)
		self.Tree.SetBranchStatus("Digi_y", 1)
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
		self.Tree.SetBranchStatus("Track_ErrorY0", 1)
		self.Tree.SetBranchStatus("Track_ErrorT0", 1)
		self.Tree.SetBranchStatus("Track_missingHitLayer", 1)
		self.Tree.SetBranchStatus("Track_expectedHitLayer", 1)
		self.Tree.SetBranchStatus("track_ipDistance", 1)
		self.Tree.SetBranchStatus("Track_hitIndices", 1)
		self.Tree.SetBranchStatus("Track_beta", 1)


	def Trigger(self):
		if len(self.Tree.Digi_x) < 3:
			return False 
		return True

	def Selection(self):
		
		if not self.Trigger():
			return False

		###########################################
		#counting total events
		self.events_passing_cuts[0] += 1.0
		self.events_passing_cuts_byfile[0] += 1.0
		###########################################
		
		###########################################
		#ntracks cut
		if (self.Tree.NumTracks < 2):
			return False

		self.events_passing_cuts[1] += 1.0
		self.events_passing_cuts_byfile[1] += 1.0
		###########################################

		###########################################
		#nvertices cut
		if self.Tree.NumVertices == 0:
			return False

		self.events_passing_cuts[2] += 1.0
		self.events_passing_cuts_byfile[2] += 1.0
		###########################################

		###########################################
		#fiducial vertex cut
		if not self.det.inBox(self.Tree.Vertex_x[0], self.Tree.Vertex_y[0], self.Tree.Vertex_z[0]):
			return False

		self.events_passing_cuts[3] += 1.0
		self.events_passing_cuts_byfile[3] += 1.0
		###########################################

		###########################################
		#floor veto w/ expected hit cuts
		for hity in self.Tree.Digi_y:
			if self.det.inLayer(hity) < 2:
				return False

		expected_hits = util.unzip(self.Tree.Track_expectedHitLayer)

		bottom_layer_expected_hits = []

		for exp_list in expected_hits:
			for val in exp_list:
				if val < 2:
					bottom_layer_expected_hits.append(val)

		if len(bottom_layer_expected_hits) < 3:
			return False


		self.events_passing_cuts[4] += 1.0
		self.events_passing_cuts_byfile[4] += 1.0
		###########################################

		###########################################
		#vertex before track cut

		
		return True


	def Plot1D(self, branch_name):
		plotvar = 0.0
		self.Tree.SetBranchStatus("")













