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

if __name__ == "__main__":

    warnings.filterwarnings("ignore")


    '''
    RUN OPTION
    '''
    opt = 13



    if opt == 0: # visualise

        bracket = [34, 47, 70, 81, 96, 147, 176, 215, 234, 
		251, 268, 283, 295, 319, 395, 398, 404, 412, 
                422, 432, 448, 487, 494, 569, 570, 581, 587, 
                599, 605, 608, 629, 633, 642, 643, 674, 678, 
                735, 898, 920, 926, 951, 990]

        for evnum in bracket:
            ev = event.Event(sys.argv[1], evnum)
            ev.ExtractTruthPhysics()
            ev.Print()

            ev.kalman = True
            ev.scattered = True
            ev.merged = False

            ev.unused = True
            ev.used = True
            ev.vertex = True
            ev.merged = True
            ev.DrawRecoPoint()

            ev.merged = True

            ev.DrawRecoPoint()

    elif opt == 1:
        ev.PlotMomentum(num=1000, cut=1000)

    elif opt == 2:
        ev.PlotHitAccuracy(num=100)

    elif opt == 3:
        ev.kalman = True

        ev.DrawRecoPoint()

    elif opt == 4:

        ev = event.Event(sys.argv[1], 0)
        missed_inds_k, missed_inds_lin, high_mult = ev.find_missed_tracks(num=(0,1000), mult=3)

        # print(missed_inds_k)
        # print("Missed by kalman: ", len(missed_inds_k))
        # print(missed_inds_lin)
        # print("Missed by lin: ", len(missed_inds_lin))

        print(high_mult)


    elif opt == 5: # for automator

        ev = event.Event(sys.argv[1], 0)
        missed_inds_k, missed_inds_lin, high_mult = ev.find_missed_tracks(num=(0,1000), mult=10)

        '''
        print(missed_inds_k)
        print("Missed by kalman: ", len(missed_inds_k))
        print(missed_inds_lin)
        print("Missed by lin: ", len(missed_inds_lin))
        '''
        print("high multiplicity events: \n",high_mult)

        for evn in high_mult:

            ev = event.Event(sys.argv[1], evn)
#            ev = event.Event(file, evn)
            ev.writeDirectory = str(sys.argv[2])

            ev.unused = True
            ev.used = True
            ev.vertex = True

            ev.ExtractTruthPhysics()
            ev.Print()

            ev.kalman = True

            ev.DrawRecoPoint()




#            ev.kalman = False

#            ev.DrawRecoPoint()

    elif opt == 6: # for automator
        # missed_inds_k, missed_inds_lin, high_mult = ev.find_missed_tracks(num=(0,5000), mult=1)
        ev = event.Event(sys.argv[1], 0)
        bracket = ev.find_bracket(num=(0,5000))

        print(bracket)
        print(len(bracket))
        # bracket = [542, 867]
        ev.writeDirectory = str(sys.argv[2])



        for evn in bracket:
            print("Current Event:", evn)

            ev = event.Event(sys.argv[1], evn)

            ev.unused = True
            ev.used = True
            ev.vertex = True

            ev.ExtractTruthPhysics()
            ev.Print()
            ev.kalman = True
            ev.merged = False

            ev.DrawRecoPoint()
            ev.merged = True


            ev.DrawRecoPoint()

    elif opt == 7:

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

    elif opt == 8:

        inds = []

        for ind in inds:
            ev = event.Event(sys.argv[1],ind)
            ev.writeDirectory = str(sys.argv[2])

            ev.ExtractTruthPhysics()
            ev.Print()


            if ev.num_truth < 2:
                ev.kalman = True
                ev.merged = True
                ev.vertex = True
                ev.DrawRecoPoint()

    elif opt == 9: # plot vertex info

        def bool_func(tree, op=0):

            tree.SetBranchStatus("Vertex_k_m_x",1)
            tree.SetBranchStatus("Vertex_x",1)

            tree.SetBranchStatus("NumTracks_k_m",1)
            tree.SetBranchStatus("NumTracks_k_m",1)

            tree.SetBranchStatus("king_move_inds",1)

            if op == 0 and len(tree.Vertex_k_m_x) >= 1: # event selection criteria
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

            ev.vert_vel = False
            ev.vertex = True
            ev.kalman = False
            ev.DrawRecoPoint()

            ev.VertexAccuracy()

    elif opt == 10:

        ev = event.Event(sys.argv[1], 0)
        ev.PlotBeta()

    elif opt == 11:
        ev = event.Event(sys.argv[1], 0)
        ev.PlotVertexBeta()

    elif opt == 12:

        ev = event.Event(sys.argv[1], 0)

        drs = []

        max = 1e3

        for ev_num in range(ev.Tree.GetEntries()):

            print("event {}".format(ev_num)) if ev_num % 100 == 0 else None

            ev = event.Event(sys.argv[1], ev_num)

            dr = ev.VertexAccuracy()

            drs.extend(dr)

        vis_engine = visualization.Histogram(drs, rng=(0,max), Title='Vertex Resolution', xaxis='Vertex Distance from Decay Location [cm]')


    elif opt == 13:

        ev = event.Event(sys.argv[1], 0)
        ev.TrackEfficiency()
