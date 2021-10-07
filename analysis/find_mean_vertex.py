#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import ROOT as root
import sys
import util
import warnings
import numpy as np
import matplotlib.pyplot as plt

# [x, y, z] in CMS

# [y, -z + 80, x] * 100 in Geant
#script_vert_loc = [0, (80 + 20) * 100, 120 * 100] # [cm]

# [y, -z + 85.47, x] * 100 in Geant
script_vert_loc = [0, (85.47 + 20) * 100, 120 * 100] # [cm]


def main():

	file = sys.argv[1]

	tfile = root.TFile.Open(file)
	Tree = tfile.Get("integral_tree")

	Tree.SetBranchStatus("Track_k_m_velX", 1)
	Tree.SetBranchStatus("Track_k_m_velY", 1)
	Tree.SetBranchStatus("Track_k_m_velZ", 1)
	Tree.SetBranchStatus("Track_k_m_x0", 1)
	Tree.SetBranchStatus("Track_k_m_y0", 1)
	Tree.SetBranchStatus("Track_k_m_z0", 1)

	Tree.SetBranchStatus("Hit_particlePx",1)
	Tree.SetBranchStatus("Hit_particlePy",1)
	Tree.SetBranchStatus("Hit_particlePz",1)
	Tree.SetBranchStatus("Hit_x",1)
	Tree.SetBranchStatus("Hit_y",1)
	Tree.SetBranchStatus("Hit_z",1)

	ys, zs = [], []
	ys_t, zs_t = [], []

	for ev in range(Tree.GetEntries()):
		Tree.GetEntry(ev)

		for tr in range(len(Tree.Track_k_m_x0)):
			# use x here since there is no issue with it
			dx = script_vert_loc[0] - Tree.Track_k_m_x0[tr]
			dt = dx / Tree.Track_k_m_velX[tr]

			y_pred = Tree.Track_k_m_y0[tr] + Tree.Track_k_m_velY[tr] * dt
			z_pred = Tree.Track_k_m_z0[tr] + Tree.Track_k_m_velZ[tr] * dt

			ys_t.append(y_pred)
			zs_t.append(z_pred)

		max_y = 0
		max_ind = -1

		for ht in range(len(Tree.Hit_x)):
			if Tree.Hit_y[ht] > max_y and Tree.Hit_particlePx[ht] != 0:
				max_y = Tree.Hit_y[ht]
				max_ind = ht

		tr = max_ind

		if tr != -1:

			# use x here since there is no issue with it
			dx = script_vert_loc[0] - Tree.Hit_x[tr]
			dt = dx / Tree.Hit_particlePx[tr]

			y_pred = Tree.Hit_y[tr] + Tree.Hit_particlePy[tr] * dt
			z_pred = Tree.Hit_z[tr] + Tree.Hit_particlePz[tr] * dt

			ys.append(y_pred)
			zs.append(z_pred)

	plot(ys, zs, False)
	plot(ys_t, zs_t)



def plot(ys, zs, track=True):

	ys = np.array(ys)
	zs = np.array(zs)

	#print(ys)
	#print(zs)

	inds = np.concatenate((np.where(zs > 14000)[0], np.where(zs < 10000)[0], np.where(ys < 8000)[0], np.where(ys > 12000)[0]))

	#print(inds)

	ys = np.delete(ys, inds)
	zs = np.delete(zs, inds)

	fig, ax = plt.subplots(figsize=(8,5))

	plt.scatter(ys, zs, c='k', s=3, label='predicted')
	plt.scatter(np.mean(ys), np.mean(zs), c='r', s=20, marker='x', label='mean predicted')
	plt.scatter(script_vert_loc[1], script_vert_loc[2], c='g', s=20, marker='x', label='script vertex')

	ax.set_title("Vertex Predicted Position vs Par-Card Position ({})".format(['sim','track'][int(track)]))
	ax.set_xlabel("y [cm]")
	ax.set_ylabel("z [cm]")

	plt.legend()
	plt.savefig('VertexComparison_{}.png'.format(['sim','track'][int(track)]))



if __name__=="__main__":

	main()
