#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import ROOT as root
import sys
import warnings
import glob


class Plot:

    def __init__(self):

        self.write_dir = sys.argv[1]

        self.filelist = [filename for filename in glob.iglob(self.write_dir+'/trees/*.root', recursive=True)]


    def Make_Plots(self, branches):

        for file in self.filelist:
             infile = root.TFile.Open( file ," READ ")
             self.tree = infile.Get("integral_tree")

             #print("Plotting in file: ", file)

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
		"Track_numHits",
		"local_chi_f",
		"local_chi_s",
		"Track_k_smooth_chi_sum",
		"NumTracks",
		"NumTracks_k",
		"NumTracks_k_m",
		"vertex_k_f_beta",
		"vertex_k_s_beta",
		"Track_k_beta",
		"Track_k_beta_err"]

    plt = Plot()

    plt.Make_Plots(branches)
