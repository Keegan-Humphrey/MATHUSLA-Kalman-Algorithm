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

'''
    for file in filelist:

        infile = root.TFile.Open( file ," READ ")
        print("Plotting in file: ", file)
        tree = infile.Get("integral_tree")

        canvas = root.TCanvas()
        canvas.cd()

        tree.Draw("Track_k_numHits")
        canvas.SaveAs(write_dir+"/plots/Track_k_numHits/"+file.split('/')[-1].split('.')[0]+"_Track_k_numHits.png")

        tree.Draw("Track_numHits")
        canvas.SaveAs(write_dir+"/plots/Track_numHits/"+file.split('/')[-1].split('.')[0]+"_Track_numHits.png")

        tree.Draw("local_chi_f")
        canvas.SaveAs(write_dir+"/plots/local_chi_f/"+file.split('/')[-1].split('.')[0]+"_local_chi_f.png")

        tree.Draw("local_chi_s")
        canvas.SaveAs(write_dir+"/plots/local_chi_s/"+file.split('/')[-1].split('.')[0]+"_local_chi_s.png")

        tree.Draw("Track_k_chi2PerNdof")
        canvas.SaveAs(write_dir+"/plots/Track_k_chi2PerNdof/"+file.split('/')[-1].split('.')[0]+"_Track_k_chi2PerNdof.png")

        tree.Draw("NumTracks")
        canvas.SaveAs(write_dir+"/plots/NumTracks/"+file.split('/')[-1].split('.')[0]+"_NumTracks.png")

        tree.Draw("NumTracks_k")
        canvas.SaveAs(write_dir+"/plots/NumTracks_k/"+file.split('/')[-1].split('.')[0]+"_NumTracks_k.png")

        tree.Draw("NumTracks_k_m")
        canvas.SaveAs(write_dir+"/plots/NumTracks_k_m/"+file.split('/')[-1].split('.')[0]+"_NumTracks_k_m.png")

        tree.Draw("vertex_k_f_beta")
        canvas.SaveAs(write_dir+"/plots/vertex_k_f_beta/"+file.split('/')[-1].split('.')[0]+"_vertex_k_f_beta.png")

        tree.Draw("vertex_k_s_beta")
        canvas.SaveAs(write_dir+"/plots/vertex_k_s_beta/"+file.split('/')[-1].split('.')[0]+"_vertex_k_s_beta.png")

        tree.Draw("Track_k_beta")
        canvas.SaveAs(write_dir+"/plots/Track_k_beta/"+file.split('/')[-1].split('.')[0]+"_Track_k_beta.png")

        tree.Draw("Track_k_beta_err")
        canvas.SaveAs(write_dir+"/plots/Track_k_beta_err/"+file.split('/')[-1].split('.')[0]+"_Track_k_beta_err.png")

        canvas.SetLogy(10)

        tree.Draw("Track_k_beta")
        canvas.SaveAs(write_dir+"/plots/Track_k_beta/"+file.split('/')[-1].split('.')[0]+"_Track_k_beta_log.png")

        tree.Draw("Track_k_beta_err")
        canvas.SaveAs(write_dir+"/plots/Track_k_beta_err/"+file.split('/')[-1].split('.')[0]+"_Track_k_beta_err_log.png")

        tree.Draw("vertex_k_f_beta")
        canvas.SaveAs(write_dir+"/plots/vertex_k_f_beta/"+file.split('/')[-1].split('.')[0]+"_vertex_k_f_beta_log.png")

        tree.Draw("vertex_k_s_beta")
        canvas.SaveAs(write_dir+"/plots/vertex_k_s_beta/"+file.split('/')[-1].split('.')[0]+"_vertex_k_s_beta_log.png")

        tree.Draw("NumTracks")
        canvas.SaveAs(write_dir+"/plots/NumTracks/"+file.split('/')[-1].split('.')[0]+"_NumTracks_log.png")

        tree.Draw("NumTracks_k")
        canvas.SaveAs(write_dir+"/plots/NumTracks_k/"+file.split('/')[-1].split('.')[0]+"_NumTracks_k_log.png")

        tree.Draw("NumTracks_k_m")
        canvas.SaveAs(write_dir+"/plots/NumTracks_k_m/"+file.split('/')[-1].split('.')[0]+"_NumTracks_k_m_log.png")

        infile.Close()
'''
