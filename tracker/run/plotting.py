#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import ROOT as root
import sys
import warnings
import glob


class Plot:

    def __init__(self):

	# output directory for a tracker run with ./run.sh
        self.write_dir = sys.argv[1]

        self.filelist = [filename for filename in glob.iglob(self.write_dir+'/trees/*.root', recursive=True)]


    def Make_Plots(self, branches):

        # make write directories for each branch
        for branch in branches:
             directory = self.write_dir+'/plots/'+branch+'/'

             #if not os.path.exists(directory):
             os.makedirs(directory)
#             print("made "+directory)

        for file in self.filelist:
             infile = root.TFile.Open( file ," READ ")
             self.tree = infile.Get("integral_tree")

             print("Plotting in file: ", file)

             self.canvas = root.TCanvas()
             self.canvas.cd()

             for branch in branches:
                 self.Draw_Save(file, branch) # linear scale plots

             self.canvas.SetLogy(10)

             for branch in branches:
                 self.Draw_Save(file, branch, True) # log scale plots

             infile.Close()


    def Draw_Save(self, file, branch, log=False):

        self.tree.Draw(branch)

        self.canvas.SaveAs(self.write_dir \
    		+ "/plots/{}/".format(branch) \
    		+ file.split('/')[-1].split('.')[0] \
                + "_{}{}.png".format(branch,['','_log'][int(log)]))



if __name__ == "__main__":
    branches = ["Track_k_numHits",
		"local_chi_f",
		"local_chi_s",
		"Track_k_smooth_chi_sum",
		"NumTracks_k_m",
		"Track_k_beta",
		"Track_k_beta_err",
		"Vertex_chi2",
		"Track_k_m_opening_angle",
		"Hit_energy",
		"Track_k_m_x_std_scat_per_m",
		"Track_k_m_z_std_scat_per_m",
		"Digi_energy"]

    '''
    branches = ["Vertex_k_m_ErrorT",
		"Vertex_k_m_ErrorX",
		"Vertex_k_m_ErrorY",
		"Vertex_k_m_ErrorZ"]

    branches = ["Track_k_filterchi",
		"Track_k_smoothchi",
		"local_chi_f",
                "local_chi_s",
		"Track_k_beta"]

    '''
    branches = []

    plt = Plot()

    plt.Make_Plots(branches)
