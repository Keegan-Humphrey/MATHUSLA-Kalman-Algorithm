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
import joblib


'''
RUN OPTION
'''
 
option = 0


print('hello viewer')

def main(opt):
    ''' opt (option) determines what process to run '''

    try:
        ev = event.Event(sys.argv[1], 0)

    except IndexError: # no argument passed (some opts don't use one)
        pass

    if opt == 0:
        # plot visualisations for events used in the analysis
        # to be put in ./vis_plots/
        
        cut = -1 # -1 to look at survivors (otherwise indexed as in flows)
        
        passed_events = joblib.load('passed_events.joblib')
        #passed_events = joblib.load('passed_events_1e5_8_5_22.joblib')

#        passed_events = joblib.load('passed_events_run6_4hits_23_2_22.joblib')
        
        total_event_cap = 45
        file_event_cap = 15

        survivors = True # False => look at events cut at cut, True => look at survivors

        if not os.path.exists('vis_plots'):
            os.makedirs('vis_plots')

        #files_handled = 0
        
        total_events_processed = 0

        for file in passed_events.keys():
        #file = '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/21_02_22/13_52_17/trees/stat_2_0.root'

        #if True:


            events_processed = 0
            #ev = event.Event(file, 0)
            try:

                if survivors:
                    inds = passed_events[file][cut].astype(int)

                else:
                    inds_before = set(passed_events[file][cut-1].astype(int))
                    inds_after = set(passed_events[file][cut].astype(int))
                
                    inds = inds_before - inds_after # show events cut

            except:
                inds = []

            if len(inds) != 0:
                #print(files_handled)
                print(file)
                #print(inds)

            #files_handled += 1
            
            #if file != '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm//08_01_22/17_50_20/trees/stat_0_0.root':
            #    continue
        
            for ind in inds:
                if events_processed > file_event_cap:
                    break
            
                print("Events Processed: ",total_events_processed)
                events_processed += 1
                total_events_processed += 1
    
                print("Event number: ",ind)
                
    #                ev = event.Event(sys.argv[1],ind)
                ev = event.Event(file,ind)
                
                #ev.writeDirectory = str(sys.argv[2])
                ev.writeDirectory = 'vis_plots/'
    
                ev.used = True
                ev.unused = True
                ev.vert_vel = True
                ev.kalman = True
    
                ev.ExtractTruthPhysics()
                ev.Print()
    
                ev.RecoVertex()            
                #ev.RecoLinVertex()
    
                
            if total_events_processed > total_event_cap:
                break
            

    if opt == 1:
        # plot vertex info from command line file and selection function

        def bool_func(tree, op=0):
            ''' event selection function '''

            tree.SetBranchStatus("Vertex_x",1)
            tree.SetBranchStatus("Vertex_k_m_x",1)
            tree.SetBranchStatus("NumTracks_k_m",1)

            # choose (or make a new) event selection criteria
            if op == 0 and len(tree.Vertex_k_m_x) >= 1:
                return True

            if op == 1 and tree.NumTracks_k_m >= 1:
                return True

            if op == 2 and tree.NumTracks_k_m >= 2 and len(tree.Vertex_k_m_x) == 0:
                return True

            if op == 3 and len(tree.Vertex_k_m_x) >= 1 and len(tree.Vertex_x) >= 1:
                return True

            else:
                return False

        
        file = sys.argv[1]
    
        ev = event.Event(file, 0)

        inds = ev.find_with_bool(bool_func, Op=0)
        inds = inds[:50]

        #inds = np.arange(11)
        
        
        if not os.path.exists('vis_plots'):
            os.makedirs('vis_plots')
       
        for ind in inds:
            
            print("Event number: ",ind)
            
#            ev = event.Event(sys.argv[1],ind)
            ev = event.Event(file,ind)
            
            ev.writeDirectory = str(sys.argv[2])
            #ev.writeDirectory = 'vis_plots/'

            ev.used = True
            ev.unused = True
            ev.vert_vel = True
            ev.kalman = True

            ev.ExtractTruthPhysics()
            ev.Print()

            ev.RecoVertex()            
            #ev.RecoLinVertex()


    elif opt == 2:
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
        ev.VertexAccuracy()

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

    elif opt == 10:
        ev.GunStudies()

    elif opt == 11:
        ev.Chi_by_ndof()

    elif opt == 12:
        ev.PlotExpectedPositions()

    elif opt == 13:
        hmu = analyzer.H_mumu_Analyzer('/home/keeganh/scratch/stat_files//27_08_21/09_59_36/trees/')

        hmu.Plot()

    elif opt == 14:
        ev.Compare_sim_and_digi()


if __name__ == "__main__":

    warnings.filterwarnings("ignore")

    main(option)
