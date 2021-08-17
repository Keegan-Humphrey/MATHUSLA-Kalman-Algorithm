#!/usr/bin/env python
import detector
import numpy as np
import ROOT as root
import sys

if __name__ == "__main__":

	file_name = sys.argv[1]
	tfile = root.TFile.Open(file_name)
	tree = tfile.Get("box_run")
	det = detector.Detector()
	gen_trigger = 0.0
	for n in range(tree.GetEntries()):
		tree.GetEntry(n)
		n_good_tracks = 0.0
		for k, index in enumerate(tree.GenParticle_G4index):
			if index > 0:
				pdg = tree.GenParticle_pdgid[k]
				
				if pdg == 111:
					continue
				if pdg == 2112:
					continue
				if pdg == 2114:
					continue

				x, y, z = tree.GenParticle_y[k]/10., tree.GenParticle_x[k]/10., tree.GenParticle_z[k]/10.
				px, py, pz = tree.GenParticle_py[k], tree.GenParticle_px[k], tree.GenParticle_pz[k]
				#print(x, y, z)
				#print(px, py, pz)
				#print(det.nSensitiveLayers(x, y, z, px, py, pz))
				if det.nSensitiveLayers(x, y, z, px, py, pz) >= 4:
					n_good_tracks += 1.0

		if n_good_tracks >= 3:
			gen_trigger += 1.0


	print("passed gen trigger: " + str(gen_trigger))




