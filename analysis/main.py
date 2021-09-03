#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import event
import ROOT as root
import analyzer
import sys
import util
import warnings
import numpy as np
import visualization


'''
RUN OPTION
'''

option = 9



def main(opt):
    ''' opt (option) determines what process to run '''

    ev = event.Event(sys.argv[1], 0)


    elif opt == 4:
        # plot vertex info

        def bool_func(tree, op=0):
            ''' event selection function '''

            tree.SetBranchStatus("Vertex_k_m_x",1)
            tree.SetBranchStatus("NumTracks_k_m",1)

            # choose (or make a new) event selection criteria
            if op == 0 and len(tree.Vertex_k_m_x) >= 1:
                return True

            if op == 1 and tree.NumTracks_k_m == 1:
                return True

            if op == 2 and tree.NumTracks_k_m >= 2 and len(tree.Vertex_k_m_x) == 0:
                return True

            else:
                return False

        ev = event.Event(sys.argv[1], 0)
        inds = ev.find_with_bool(bool_func, Op=0)

        print(inds)

        for ind in inds:

            print("Event number: ",ind)
            ev = event.Event(sys.argv[1],ind)
            ev.writeDirectory = str(sys.argv[2])

            ev.used = True
            ev.unused = True
            ev.vert_vel = True
            ev.kalman = True

            ev.ExtractTruthPhysics()
            ev.Print()

            ev.RecoVertex()


    elif opt == 3:
        ''' find number of truth tracks per event '''

        num_truths = []

        for i in range(5000):
            ev = event.Event(sys.argv[1], i)

            ev.writeDirectory = str(sys.argv[2])

            ev.ExtractTruthPhysics()
            ev.Print()

            num_truths.append(ev.num_truth)

        num_truths = np.array(num_truths)

        num_truths[num_truths < 2] = 0
        num_truths = np.count_nonzero(num_truths)

        print("number of truth tracks \n",num_truths)


    elif opt == 3:

#        drs = []
        dx, dy, dz = [], [], []

        max = 1e3

        for ev_num in range(ev.Tree.GetEntries()):

            print("event {}".format(ev_num)) if ev_num % 100 == 0 else None

            ev = event.Event(sys.argv[1], ev_num)

#            dr = ev.VertexAccuracy()
            _dx, _dy, _dz = ev.VertexAccuracy()

             #drs.extend(dr)
            dx.extend(_dx)
            dy.extend(_dy)
            dz.extend(_dz)

#        vis_engine = visualization.Histogram(drs, rng=(0,max), Title='Vertex Resolution', xaxis='Vertex Distance from Decay Location [cm]')
        vis_engine = visualization.Histogram(dx, rng=(-max,max), Title='Vertex X Resolution', \
		xaxis='Vertex X - Decay Location X [cm]', fname='resolutionX.png')
        vis_engine = visualization.Histogram(dy, rng=(-max,max), Title='Vertex Y Resolution', \
		xaxis='Vertex Y - Decay Location Y [cm]', fname='resolutionY.png')
        vis_engine = visualization.Histogram(dz, rng=(-max,max), Title='Vertex Z Resolution', \
		xaxis='Vertex Z - Decay Location Z [cm]', fname='resolutionZ.png')

    elif opt == 4:
        ev.PlotMomentum(num=1000, cut=1000)

    elif opt == 6:
        ev.PlotBeta()

    elif opt == 7:
        ev.PlotVertexBeta()

    elif opt == 8:
        ev.TrackEfficiency()

    elif opt == 9:
        ev.TrackResolution()



if __name__ == "__main__":

    warnings.filterwarnings("ignore")

    main(option)
