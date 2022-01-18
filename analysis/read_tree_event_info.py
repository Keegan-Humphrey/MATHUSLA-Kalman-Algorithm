#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import physics
import visualization
from detector import Detector
import util
import ROOT as root


class Read:

    def __init__(self) 

#        self.find_dir = files

#        self.filelist = [filename for filename in \
#		            glob.iglob(self.find_dir+'/*.root', recursive=True)]

        self.files = 
        self.events =


    def Show_Branches(self, events, branches):

        for file in self.filelist:
             infile = root.TFile.Open( file ," READ ")

#             self.tree = infile.Get("box_run")
             self.tree = infile.Get("integral_tree")

             self.canvas = root.TCanvas()
             self.canvas.cd()

             for branch in branches:
                 self.Gather(events, branch)

             infile.Close()

             self.Show()


    def Gather(self, file, branch, log=False):

        self.tree.Draw(branch)

        self.canvas.SaveAs(self.write_dir \
                + "{}_{}.png".format(file.split('/')[-1].split('.')[0], branch))


    def Show(self):
    
    
    
if __name__=="__main__":

    branches = ['Vertex_k_m_t',
                'Digi_time',
                'Digi_y']

    reader = Read()
    
    reader.Show_Branches(branches)






    
    
    