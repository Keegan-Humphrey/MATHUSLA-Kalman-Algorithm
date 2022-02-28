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

        self.filelist = [filename for filename in \
		glob.iglob(self.write_dir+'/*.root', recursive=True)]


    def Make_Plots(self, branches):

        for file in self.filelist:
             infile = root.TFile.Open( file ," READ ")

#             self.tree = infile.Get("box_run")
             self.tree = infile.Get("integral_tree")

             self.canvas = root.TCanvas()
             self.canvas.cd()

             for branch in branches:
                 self.Draw_Save(file, branch) # linear scale plots

             infile.Close()


    def Draw_Save(self, file, branch, log=False):

        self.tree.Draw(branch)

        self.canvas.SaveAs(self.write_dir \
                + "{}_{}.png".format(file.split('/')[-1].split('.')[0], branch))



if __name__ == "__main__":

    '''
    branches = ["x_estimates_m:z_estimates_m", \
		"y_estimates_m:z_estimates_m", \
		"x_estimates_m:y_estimates_m"]

    '''
#    branches = ["Hit_x:Hit_z"] #, \
#		"Hit_y:Hit_z", \
#		"Hit_x:Hit_y"]
    '''
    branches = ["Track_k_x0:Track_k_z0", \
		"Track_k_y0:Track_k_z0", \
		"Track_k_x0:Track_k_y0"]
    '''
#    branches_1 = ["Digi_x:Digi_z"] #, \
#                "Digi_y:Digi_z", \
#                "Digi_x:Digi_y"]

#    branches.extend(branches_1)

    branches = ['Hit_x',
		'Hit_y',
		'Hit_z']

    plt = Plot()

    plt.Make_Plots(branches)
