
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import matplotlib.pyplot as plt
from itertools import combinations
from scipy import constants
from math import acos
import visualization
import physics
import ROOT as root
import numpy as np
import detector
import util
import time
import scipy.optimize as spo
import sys
import glob
import pandas as pd
import joblib
import inspect    

ncuts = 14



class sample_space():


    def __init__(self,file):

        self.file = file
        self.det = detector.Detector()
        self.cuts = dict() # information about how to use the cuts

    def inside_box(self, x, y, z):

        box_lims = self.det.BoxLimits
        y_floor = self.det.LayerYLims[2][1]

        if box_lims[0][0] < x and box_lims[0][1] > x:
           if y_floor < y and box_lims[1][1] > y: # only accept it if above the floor layers
                if box_lims[2][0] < z and box_lims[2][1] > z:
                    return True

        return False


    def in_layer(self,y_val):

        LAYERS_Y = self.det.LayerYLims
        for n in range(len(LAYERS_Y)):
            _min = LAYERS_Y[n][0]
            _max = LAYERS_Y[n][1]
            if (y_val > _min) and (y_val < _max):
                    return n

        return 999


    def get_states_m(self):

        global ncuts

        tracking_file = root.TFile.Open(self.file)
        self.tree = tracking_file.Get("integral_tree")
        self.tree.SetBranchStatus("*", 0)

        self.tree.SetBranchStatus("Digi_x", 1)
        self.tree.SetBranchStatus("Digi_y", 1)
        self.tree.SetBranchStatus("Digi_z", 1)
        self.tree.SetBranchStatus("Digi_time", 1)
        self.tree.SetBranchStatus("NumTracks_k_m", 1)

        self.tree.SetBranchStatus("Vertex_k_m_t", 1)
        self.tree.SetBranchStatus("Vertex_k_m_x", 1)
        self.tree.SetBranchStatus("Vertex_k_m_y", 1)
        self.tree.SetBranchStatus("Vertex_k_m_z", 1)
        self.tree.SetBranchStatus("Vertex_k_m_ErrorT", 1)
        self.tree.SetBranchStatus("Vertex_k_m_ErrorX", 1)
        self.tree.SetBranchStatus("Vertex_k_m_ErrorY", 1)
        self.tree.SetBranchStatus("Vertex_k_m_ErrorZ", 1)
        self.tree.SetBranchStatus("Vertex_k_m_trackIndices", 1)
        self.tree.SetBranchStatus("NumVertices_k_m", 1)

        self.tree.SetBranchStatus("Track_k_m_x0", 1)
        self.tree.SetBranchStatus("Track_k_m_velX", 1)
        self.tree.SetBranchStatus("Track_k_m_velY", 1)
        self.tree.SetBranchStatus("Track_k_m_velZ", 1)
        self.tree.SetBranchStatus("Track_k_m_hitIndices", 1)
        self.tree.SetBranchStatus("Track_k_m_expected_hit_layer", 1)

        self.tree.SetBranchStatus("Track_k_m_x0", 1)
        self.tree.SetBranchStatus("Track_k_m_y0", 1)
        self.tree.SetBranchStatus("Track_k_m_z0", 1)
        self.tree.SetBranchStatus("Track_k_smooth_chi_sum",1)


        cut_vectors = []
        
        self.additional_info = np.zeros((int(self.tree.GetEntries()),2)) # Fill me up with info if you want internal 
                                                                         # access to the cuts to plot something new
        for event_number in range(int(self.tree.GetEntries())):

            self.event_number = event_number

            self.current_vector = np.zeros(1+ncuts) # first will be the event number
            self.n = 0

            self.tree.GetEntry(event_number)
            if event_number % 1000 == 0:
                print("event:", event_number)

            #------------- record the event number for reference
            self.current_vector[0] = event_number
            self.n += 1

            #------------- cut on the number of made tracks in the event
            
            self.Tracks()
            self.n += 1

            #------------- cut on the number of vertices made in the event
            
            self.Vertices()
            self.n += 1
            
            #------------- Optional Cuts
            
            self.Fiducial_vertex()
            self.n += 1
            
            self.Floor_hits_before_vertex()
            self.n += 1

            self.Expected_hits_old()
            self.n += 1 
            
            self.Opening_angle()
            self.n += 1

            self.Topological()
            self.n += 1

            self.Vertex_track_beta()
            self.n += 1

            self.Track_hit_diffs()
            self.n += 1

            self.Exp_hits()
            self.n += 1

            self.No_floor_hits()
            self.n += 1
            
            self.Missing_hit_sum()
            self.n += 1
            
            self.Chi_ndof_cut()
            self.n += 1

            #------------ All done
            cut_vectors.append(self.current_vector)

        self.cut_vectors = np.array(cut_vectors)

        return self.cut_vectors


    def Tracks(self):
        self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'} # record the index of the cut in cut vector
                                                                            # and how to operate on the value.
                                                                            # Key is the name of this function
        self.current_vector[self.n] = self.tree.NumTracks_k_m
            
            
        
    def Vertices(self):
        self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'} 
        
        self.current_vector[self.n] = len(self.tree.Vertex_k_m_x)       
        

    def Fiducial_vertex(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'} 
            
            #------------- Fiducial vertex cut
            inside = False

            for num in range(len(self.tree.Vertex_k_m_x)):
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]
                if self.inside_box(vtxx, vtxy, vtxz):
                    inside = True

            # check if at least one vertex is in the detector
            self.current_vector[3] = int(inside) # fiducial
            
            

    def Floor_hits_before_vertex(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'} 
           
            #------------- Floor hits before vertex
            floorveto = False

            vertex_trackIndices = self.tree.Vertex_k_m_trackIndices
            self.vertex_trackIndices = util.unzip(vertex_trackIndices)

            min_by_verts = []

            for num in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vertexveto = False
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]

                ''' x-z distance of hit '''
                dists = []
                ''' ----------- '''

                for hit in range(len(self.tree.Digi_y)): # loop over the hits
                    if self.in_layer(self.tree.Digi_y[hit]) <= 2 or self.tree.Digi_z[hit] < self.det.ModuleZLims[0][0]: # check if layer index is a floor layer index or in the wall
                        if self.tree.Vertex_k_m_t[num] + self.tree.Vertex_k_m_ErrorT[num] > self.tree.Digi_time[hit]: # check if the floor or wall hit happened before the vertex
                            vertexveto = True
                            
                            ''' //----------------------------// '''
                            expected_x = []
                            expected_z = []
            
                            for ind in self.vertex_trackIndices[num]: # track
                                    ind = int(ind)
            
                                    x = self.tree.Track_k_m_x0[ind]
                                    y = self.tree.Track_k_m_y0[ind]
                                    z = self.tree.Track_k_m_z0[ind]
            
                                    vx = self.tree.Track_k_m_velX[num]
                                    vy = self.tree.Track_k_m_velY[num]
                                    vz = self.tree.Track_k_m_velZ[num]
            
                                    del_t = (self.tree.Digi_y[hit] - y) / vy
            
                                    expected_x.append(x + vx * del_t)
                                    expected_z.append(z + vz * del_t)
                                    
                            expected_x = np.array(expected_x)
                            expected_z = np.array(expected_z)
                
                            if len(expected_x) != 0:
                                min_dist = np.amin(np.sqrt((expected_x - self.tree.Digi_x[hit])**2 + (expected_z - self.tree.Digi_z[hit])**2)) # min dist to track for the hit
                                dists.append(min_dist)
                
                if len(dists) != 0:
                    min_by_verts.append(np.amin(dists)) # min dist between a track and hit in the event that happened before the vertex
                    
                ''' //-----------------------------------// '''

                if not self.inside_box(vtxx, vtxy, vtxz): # make sure the vertex is fiducial (above the floor)
                    vertexveto = True

                if not vertexveto: # There's a vertex in the event that isn't vetoed
                    break

            if len(self.tree.Vertex_k_m_x) != 0:
                if vertexveto:
                    floorveto = True

            self.current_vector[self.n] = int(not floorveto) # cut if floorveto == True => 0 == current_vector[4] (< 1 is true)

            
            ''' \\-----------------------------------\\ '''
            if len(min_by_verts) != 0:
                self.additional_info[self.event_number][0] = np.amin(min_by_verts)
                ''' \\-----------------------------------\\ '''
                

    def Expected_hits_old(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
            
            #---------- Expected Hits (old) - deprecated / doesn't run
            '''
            expectedveto = False

            fullhitindices = util.unzip(self.tree.Track_k_m_hitIndices)
            fullexp_layers = util.unzip(self.tree.Track_k_m_expected_hit_layer)

            for j in range(int(self.tree.NumVertices_k_m)):
                bottomlayer_exp = []
                bottomlayer_hits = []
                trackindices = self.fulltrackindices[j]

                for track in trackindices:
                    exp_layers = fullexp_layers[int(track)]  # layers the track is expected to hit

                    for hit in exp_layers:
                        if hit <= 2:
                            bottomlayer_exp.append(hit) # a floor layer hit is expected, record the layer

                    hitindices = fullhitindices[int(track)]

                    for i in hitindices:
                        hity = self.tree.Digi_y[i]

                        if self.in_layer(hity) <= 2:
                                bottomlayer_hits.append(hity) # hit is in a floor layer

                # we get 3 or more expected hits, but no hits, veto
                if (len(bottomlayer_hits) < 1 and len(bottomlayer_exp) >= 1):
                    expectedveto = True
            
            self.current_vector[self.n] = int(expectedveto) # bottom layer hits
            '''
                  


    def Opening_angle(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
    
            #------------- Track Opening Angle Cut
            open_angles = []
            self.fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            
            for k1 in range(len(self.tree.Vertex_k_m_x)):
                trackindices = self.fulltrackindices[k1]
                combolist = list(combinations(trackindices, 2))

                for combo in combolist:
                    combo = [int(combo[0]),int(combo[1])]

                    cos_opening_angle = self.tree.Track_k_m_velX[combo[0]]*self.tree.Track_k_m_velX[combo[1]] + \
                                        self.tree.Track_k_m_velY[combo[0]]*self.tree.Track_k_m_velY[combo[1]] + \
                                        self.tree.Track_k_m_velZ[combo[0]]*self.tree.Track_k_m_velZ[combo[1]]
                    mag1 = np.sqrt(self.tree.Track_k_m_velX[combo[0]]**2 + self.tree.Track_k_m_velY[combo[0]]**2 + self.tree.Track_k_m_velZ[combo[0]]**2)
                    mag2 = np.sqrt(self.tree.Track_k_m_velX[combo[1]]**2 + self.tree.Track_k_m_velY[combo[1]]**2 + self.tree.Track_k_m_velZ[combo[1]]**2)
                    cos_opening_angle = cos_opening_angle/( mag1*mag2 )
                    open_angles.append(acos(cos_opening_angle))

            # sort in terms of descending opening angle
            open_angles.sort(reverse = True)
            
            # sort in terms of ascending opening angle
#            open_angles.sort(reverse = False)

            if len(open_angles) != 0:
                # cut on largest opening angle
                self.current_vector[self.n] = open_angles[0]
                  
                                                  

    def Topological(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}

            #---------------- check if vertices are below track hits within some sigma (Topological Cut)
            topoveto = 0

            self.fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            self.fullhitindices = util.unzip(self.tree.Track_k_m_hitIndices)

            pulls = []

            for k1 in range(len(self.tree.Vertex_k_m_x)):
                vertexveto = False
                trackindices = self.fulltrackindices[k1]

                min_y = 1e6
                for t1 in range(len(self.tree.Track_k_m_x0)):
                    hitindices = self.fullhitindices[t1]
                    for hit in hitindices:
                         if self.tree.Digi_y[hit] < min_y:
                             min_y = self.tree.Digi_y[hit]

                if self.tree.Vertex_k_m_ErrorY[k1] != 0:
                    pulls.append((min_y - self.tree.Vertex_k_m_y[k1]) / self.tree.Vertex_k_m_ErrorY[k1]) # pull of lowest hit from vertex
                else:
                    pulls.append(1e6)

            if len(pulls) != 0:
                self.current_vector[self.n] = np.max(pulls)

            

    def Vertex_track_beta(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
             
            #---------- Vertex Beta Cut
            # want at least one vertex with two betas within some interval in each event

            second_closest_betas = []

            for k1 in range(len(self.tree.Vertex_k_m_x)):
                betas = []

                trackindices = self.fulltrackindices[k1]

                for ind in trackindices:
                    ind = int(ind)

                    beta = np.linalg.norm([self.tree.Track_k_m_velX[ind],
		              			self.tree.Track_k_m_velY[ind],
				              	self.tree.Track_k_m_velZ[ind]]) / physics.c

                    betas.append(np.abs(1-beta)) # residue from beta = 1

                betas.sort(reverse = False) # put in ascending order

                second_closest_betas.append(betas[1]) # choose second closest

            second_closest_betas.sort(reverse = False) # put in ascending order

            if len(second_closest_betas) != 0:
#                current_vector[8] = 1 / second_closest_betas[0] 
                self.current_vector[self.n] = 1 / second_closest_betas[0] # set to vertex with closest second_closest beta
                                                                # use reciprocal so we can cut events 1/beta < 1/x (< fixed by algorithm)
                                                                

    def Track_hit_diffs(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
             
            #---------- Track Length Difference Cut
            # cut on minimum difference in number of hits in tracks in a vertex

            min_diffs = []

            for k1 in range(len(self.tree.Vertex_k_m_x)):
                diffs = []

                trackindices = self.fulltrackindices[k1]
                combolist = list(combinations(trackindices, 2))

                for combo in combolist:
                    combo = [int(combo[0]),int(combo[1])]

                    diff = np.abs(len(self.fullhitindices[combo[1]]) \
                                  - len(self.fullhitindices[combo[0]]))

                    diffs.append(diff)

                min_diffs.append(np.amin(diffs))

            if len(min_diffs) != 0:
                self.current_vector[self.n] = 1/np.amin(min_diffs) if np.amin(min_diffs) != 0 else 1e6


    def Exp_hits(self):
            self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
             
            #------------- Expected Hits Cut

            max_centre_dist = 0

            centre = [np.mean(self.det.BoxLimits[0]), np.mean(self.det.BoxLimits[2])]

            floor_ys = [np.sum(self.det.LayerYLims[i]) / 2 for i in range(3)]
            floor_y = np.sum(self.det.LayerYLims[2]) / 2

            expected_x = []
            expected_z = []

#        floor_y = np.sum(Detector.LayerYLims[2]) / 2

            hits_in_floor = False

            for i in range(3):
                    
                if np.any(self.tree.Digi_y == np.sum(self.det.LayerYLims[i]) / 2):
                    hits_in_floor = True
                    
                    #print("floor",np.sum(self.det.LayerYLims[i]) / 2)
                    #print("y hits ",np.array(self.tree.Digi_y))

            if len(self.vertex_trackIndices) > 0 and not hits_in_floor:
                for vert in range(len(self.vertex_trackIndices)): # vertex
                    vert = int(vert)

                    for ind in self.vertex_trackIndices[vert]: # track
                        ind = int(ind)

                        x = self.tree.Track_k_m_x0[ind]
                        y = self.tree.Track_k_m_y0[ind]
                        z = self.tree.Track_k_m_z0[ind]

                        vx = self.tree.Track_k_m_velX[ind]
                        vy = self.tree.Track_k_m_velY[ind]
                        vz = self.tree.Track_k_m_velZ[ind]

                        del_t = (floor_y - y) / vy

                        expected_x.append(x + vx * del_t)
                        expected_z.append(z + vz * del_t)

            edge_dist_x = (self.det.BoxLimits[0][1] - self.det.BoxLimits[0][0]) / 2 - np.abs(np.array(expected_x) - centre[0]) # distance from expected hits
            edge_dist_z = (self.det.BoxLimits[2][1] - self.det.BoxLimits[2][0]) / 2 - np.abs(np.array(expected_z) - centre[1]) # to edge of detector

            #print('dists are')
            #print(edge_dist_x)
            #print(edge_dist_z)

            if len(edge_dist_x) != 0:
                 #min_edge_dist = np.amin(np.amin([edge_dist_x, edge_dist_z], axis=1))
                 #self.current_vector[self.n] = min_edge_dist

                 # this is choosing the most negative values, no good! we want the least negative values
                 
                 #print(min_edge_dist)

                 #self.additional_info_1[self.event_number] = expected_x[:2] #edge_dist_x[:2]
                 #self.additional_info_2[self.event_number] = expected_z[:2] #edge_dist_z[:2]

                 #max_edge_dist = np.amin(np.max([edge_dist_x, edge_dist_z], axis=1))
                 #self.current_vector[self.n] = max_edge_dist
                 
                 average_edge_dist = np.amin(np.mean([edge_dist_x, edge_dist_z],axis=1)) # calculate average position, find distance to the edge
                 self.current_vector[self.n] = average_edge_dist
                 #print(max_edge_dist)
                 

    def No_floor_hits(self):
        self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
             
        ''' Veto the event if it has floor hits '''
    
        floor_veto = False
    
        for hit in range(len(self.tree.Digi_y)): # loop over the hits
            if self.in_layer(self.tree.Digi_y[hit]) <= 2: # check if layer index is a floor layer index
                floor_veto = True
                break

        if floor_veto:
            self.current_vector[self.n] = True  
            
            
    def Missing_hit_sum(self):
        self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
             
        ''' Veto the event if there are too many missed layers for the tracks making up the vertex => likely fake'''
        
        missed_hits_event = []
        
        for vert in range(len(self.vertex_trackIndices)): # vertex
            vert = int(vert)

            missed_hits_vertex = []

            expected_layers = util.unzip(self.tree.Track_k_m_expected_hit_layer)

            for trk in self.vertex_trackIndices[vert]: # track
                trk = int(trk)
    
                y = self.tree.Track_k_m_y0[trk]
    
                lowest_hit_layer = self.in_layer(y)      
    
                try:
                    expected_layers_above = len([layer for layer in expected_layers[trk] if layer >= lowest_hit_layer])
        
                except IndexError:
                    print(self.vertex_trackIndices[vert])
                    print()
                    print(expected_layers)
                    print(len(expected_layers))
                    print(trk)
        
                num_hits = len(util.unzip(self.tree.Track_k_m_hitIndices)[trk])
                
                missed_layers =  expected_layers_above - num_hits

                missed_hits_vertex.append(missed_layers)
                
            missed_hits_vertex.sort(reverse=True) # sort in ascending order of missed hits
            
            missed_hits_sum = np.sum(missed_hits_vertex[:2]) # add missed hits in tracks with the fewest

            missed_hits_event.append(missed_hits_sum)
            
        if len(missed_hits_event) != 0:
            self.current_vector[self.n] = -np.max(missed_hits_event) # if missed hit sum of the event is too large cut the event            


    def Chi_ndof_cut(self):
        self.cuts[inspect.stack()[0][3]] = {'index':self.n, 'cut if':'<'}
    
        chi_ndofs = []
        
        for vert in range(len(self.vertex_trackIndices)): # vertex
            vert = int(vert)
            
            vert_chi = []
            
            for trk in self.vertex_trackIndices[vert]: # track
                trk = int(trk)
                
                vert_chi.append(self.tree.Track_k_smooth_chi_sum[trk])
                
            vert_chi.sort(reverse=False)
            
            chi_ndofs.append(np.sum(vert_chi[:2])) # add the two lowest chi ndof in the vertex
            
        if len(chi_ndofs) != 0:
            self.current_vector[self.n] = -np.amin(chi_ndofs) # choose the lowest chis
            self.n += 1
            
            self.current_vector[self.n] = -np.amin(chi_ndofs) # choose the lowest chis
        


class scissors():
    '''since this is what does the cutting!'''

    def __init__(self,file):

        self.file = file

    def store_space(self):

        space = sample_space(self.file)
        
        self.cut_vectors = space.get_states_m()
        self.func_dicts = space.cuts
      
        self.additional_info = space.additional_info
      
        self.events = len(self.cut_vectors)
        
        print('number of events is: ',len(self.cut_vectors))

        
    def cut_dict(self, cut_options, permutation):

        global ncuts

        self.flows = np.zeros(ncuts)
        local_vectors = np.copy(self.cut_vectors)
        self.survivor_inds = []

        if local_vectors.size != 0:
            
            for i in permutation:
            
                opt_dict = cut_options[str(i)]
                func = opt_dict['func_name']
                func_dict = self.func_dicts[func]
                
                n = func_dict['index'] # cut vector index for this cut
                
                if opt_dict['on?'] == 1: 
                    cut_parameter = opt_dict['cut parameter']
                
                    if func_dict['cut if'] == '<': # only < for now, as this is what was used previously
                                                   # add an elif below to add >, ==, etc options
                        inds = np.where(local_vectors[:,n] < cut_parameter) 
                        
                else:
                    inds = [] # don't cut any events
                     
                     
                local_vectors = np.delete(local_vectors, inds, 0)

                self.flows[i] += np.shape(local_vectors)[0]
                self.survivor_inds.append(local_vectors[:,0])

        self.survivors = self.flows[-1]
        
        

class par_func():
    ''' handle multiple files to cut and create a function that takes cut parameters and returns signal to background ratios'''

    def __init__(self,backgrounds,signals):

        if type(backgrounds) == list and type(signals) == list:
            self.b_files = backgrounds
            self.s_files = signals

        else:
            self.b_files = [backgrounds]
            self.s_files = [signals]

    def get_states(self):

        self.drawers = [[],[]]

        i = 0

        for files in [self.s_files, self.b_files]:
            for fl in files:

                print(fl)

                scissor = scissors(fl)

                scissor.store_space()
                self.drawers[i].append(scissor)

            i += 1


#    def cut(self, p0, p1, p2, p3, p4, p5):
#    def cut(self, p0, p1, p4):
    def cut(self, p4):
        ''' calculate and return signal to background ratio for all handled files'''

        i = 0
        self.scores = [[],[]]

        for drawer in self.drawers:
            for scissor in drawer:

#                scissor.cut(p0, p1, p2, p3, p4, p5)
                scissor.cut(2, 1, 1, 1, -1, p4, 1)
                self.scores[i].append(scissor.survivors / \
				      np.shape(scissor.cut_vectors)[0])

            i += 1

        return np.sum(self.scores[0]) if np.sum(self.scores[1]) == 0 \
                                      else np.sum(self.scores[0]) / np.sum(self.scores[1])


#    def get_flows(self, p0, p1, p2, p3, p4, p5):
#    def get_flows(self, p0, p1, p4):
    def get_flows(self, p4):

       flows = [[],[]]

       i = 0

       for drawer in self.drawers:
            for scissor in drawer:

#                scissor.cut(p0, p1, p2, p3, p4, p5)
#                scissor.cut(p0, p1, 1, 1, 1, p4, 1)
                scissor.cut(2, 1, 1, 1, -1, p4, 1)
                flows[i].append(scissor.flows)

            i += 1

       return flows



def optimize():
    ''' deprecated code to optimize cut parameters '''

    directory = sys.argv[1]
    b_files = [filename for filename in glob.iglob(directory+'stat_0_*.root', recursive=True)]
#    s_files = [filename for filename in glob.iglob(directory+'stat_2_*.root', recursive=True)]

    b_files = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/13_10_21/13_20_32/trees/stat_4_0.root']
    s_files = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/13_10_21/08_16_40/trees/stat_0_0.root',
            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/13_10_21/08_16_40/trees/stat_2_0.root']

    print(b_files)
    print(s_files)

    for i in range(len(b_files)):

        func = par_func(b_files[i],s_files[i])

        func.get_states()

#       bounds = [(-0.01,5.01),(-0.01,1.01),(-0.01,1.01),(-0.01,1.01),(-1.01,1.01),(-0.01,1.01)]
        bounds = [(-2,2)]

        opts = ({"minimize_every_iter":True,
                 "disp":False})

        # need the lambda function since shgo doesn't like methods!
        x = spo.shgo(lambda X: -func.cut(*X), bounds, iters=10, options=opts)

        print("pars are ",x)
        print("optimal is ",func.cut(*x['x']))

        flows = func.get_flows(*x['x'])
        print("signal flows are \n ",flows[0],"\nbackground flows are \n",flows[1])
        '''
        flows = func.get_flows(2,1,1)
        print("optimal is ",func.cut(2,1,1))
        print("signal flows are \n ",flows[0],"\nbackground flows are \n",flows[1])
        '''





def main():

    
    if len(sys.argv) > 1:
        # default set for job submission with files and write locations passed as arguments
        
        if len(sys.argv) != 4:
            print('Usage: python cutflow.py <tracker tree> <save directory> <integer identifier>')
            return
    
        plot_cut = False
        plot_add_info = False
        sum_flows = False # True <=> Background
        load = False
        save = True
        TwoD = False
        
        files = [sys.argv[1]]

    else:
        plot_cut = False
        plot_add_info = False
        sum_flows = False # True <=> Background
        load = False
        save = False
        TwoD = False # whether additional info is 2D
        
        if sum_flows: # Background
            directory_3 = '/home/keeganh/scratch/job_test/W_sample_dir/run3/18_12_21/11_24_04/trees/'
            #directory_4 = '/home/keeganh/scratch/job_test/W_sample_dir/run4/09_01_22/09_11_57/trees/'
            
            files = [filename for filename in glob.iglob(directory_3+'stat_*.root', recursive=True)]
            #files.extend([filename for filename in glob.iglob(directory_4+'stat_*.root', recursive=True)])
    
        else: # signal
            files = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_0_0.root',
                '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_2_0.root',
                '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm//08_01_22/17_50_20/trees/stat_0_0.root'][-1:]
    #            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_4_0.root']
    

    cut = 5 # which cut to plot
    info_cut = 3 # which cut we take additional info from 
    
    sum_values = []
    sum_data_x = []
    sum_data_z = []
    drawers = [] # to hold the scissors

    flows = dict()
    passed_events = dict()
    
    option = ['cut name', 'cut parameter', 'units', 'on?'] # [name, value of the cut, units of what is being cut on, whether the cut is used (1 => yes)]
    func_name = 'func_name'

    cut_options = {'0' :{option[0]:'2 Tracks'                     ,option[1]:2       ,option[2]:'tracks'      ,option[3]:1 , func_name:'Tracks' },
                   '1' :{option[0]:'Vertices'                     ,option[1]:1       ,option[2]:'verts'       ,option[3]:1 , func_name:'Vertices' },
                   '2' :{option[0]:'Fiducial Vertex'              ,option[1]:1       ,option[2]:'bool'        ,option[3]:1 , func_name:'Fiducial_vertex' },
                   '3' :{option[0]:'Floor/Wall Hits Before Vertex',option[1]:1       ,option[2]:'bool'        ,option[3]:1 , func_name:'Floor_hits_before_vertex' },
                   '4' :{option[0]:'Old Expected Hits'            ,option[1]:-1      ,option[2]:'bool'        ,option[3]:0 , func_name:'Expected_hits_old' },
                   '5' :{option[0]:'Vertex opening angle'         ,option[1]:0.02    ,option[2]:'rad'         ,option[3]:1 , func_name:'Opening_angle' },
                   '6' :{option[0]:'Topological Veto'             ,option[1]:-1e2    ,option[2]:'sigma'       ,option[3]:0 , func_name:'Topological' },
                   '7' :{option[0]:'2 Good Betas in a Vertex'     ,option[1]: 1 / 0.2,option[2]:'1 / beta res',option[3]:1 , func_name:'Vertex_track_beta' },
                   '8' :{option[0]:'Hit Differences'              ,option[1]: 1 / -1 ,option[2]:'1 / hits'    ,option[3]:0 , func_name:'Track_hit_diffs' },
                   '9' :{option[0]:'Expected hit edge distance'   ,option[1]:1000    ,option[2]:'cm'          ,option[3]:1 , func_name:'Exp_hits' },
                   '10':{option[0]:'No Floor Hits'                ,option[1]:1       ,option[2]:'bool'        ,option[3]:0 , func_name:'No_floor_hits' },
                   '11':{option[0]:'Missing Hit Sum'              ,option[1]:-6      ,option[2]:'-missed hits',option[3]:1 , func_name:'Missing_hit_sum' },
                   '12':{option[0]:'Chi sum cut 1'                ,option[1]:-20     ,option[2]:'-chi ndof'   ,option[3]:1 , func_name:'Chi_ndof_cut' },
                   '13':{option[0]:'Chi sum cut 2'                ,option[1]:-10     ,option[2]:'-chi ndof'   ,option[3]:0 , func_name:'Chi_ndof_cut' }} # func_name is the corresponding data collection function in get_states_m
                       
    for opt in option:
        flows[opt] = np.array([])

    for cut_ind in cut_options.keys():
        for i in range(len(option)):
            flows[option[i]] = np.append(flows[option[i]], cut_options[cut_ind][option[i]])
            
    flows[option[3]] = flows[option[3]].astype(int)

    permutation = [0,1,2,3,9,11,4,5,6,7,8,10,12,13] # order in which the cuts are performed
                                                    # describes a permutation of 
                                                    # (0, ..., ncuts-1); keys of cut_options

    for key in flows.keys():
        flows[key] = flows[key].take(permutation,0) # permute each entry in flows dictionary
 
    total_events = 0

    i = 0

    for file in files:

        if sum_flows or len(sys.argv) > 1:
            sample = 'W'

        else:
            sample = ['h10','h2','qq'][i] # particle type for each sample in files 

        print(file)

        scissor = scissors(file)

        if not load:
            try:
                scissor.store_space() # read the trees and gather data for cutting from each event
            
            except ReferenceError: # tree is empty
                print("Sorry, that tree doesn't have anything in it")
                i += 1
                continue
            
        if not load and np.shape(scissor.cut_vectors)[0] != 0:
            if len(sys.argv) > 1:
                save_files_dir = sys.argv[2]
                job_num = sys.argv[3]
                
                joblib.dump(scissor,save_files_dir+'/scissor_{}_{}.joblib'.format(sample,job_num)) # save the scissor objects
                                                                                                       # => only need to read data once
            elif save:
                if not os.path.exists('save_files'):
                    os.makedirs('save_files')
                    
                joblib.dump(scissor,'save_files/scissor_{}_{}.joblib'.format(sample,i)) # save the scissor objects
                                                                                                       # => only need to read data once
            
        if load:
            try:
                scissor = joblib.load('save_files/scissor_{}_{}.joblib'.format(sample,i)) # load the scissor objects so we don't
                                                                                          # need to read and process the data again
                                                                                          # (just need to cut with chosen parameters)
            except FileNotFoundError:
                print("Sorry, I couldn't find that file")
                i += 1
                continue
                
        scissor.cut_dict(cut_options, permutation)
        
        passed_events[file] = scissor.survivor_inds
        
        file_flows = scissor.flows.astype(int)
        file_flows = file_flows.take(permutation,0) # permute the flows to the order the cuts were performed
        
        if sum_flows:
            try:
                flows['Flows'] += file_flows
                             
            except:
                flows['Flows'] += file_flows

            total_events += scissor.events

        else:
            flows[sample] = file_flows
            flows[sample+' fraction'] = np.round(file_flows / scissor.events, 4)



        if plot_cut and len(scissor.survivor_inds) != 0:
            #------- Plot the Chosen cut
            inds = np.array(scissor.survivor_inds[cut-1],dtype=int) # surviving indices before cut

#            values = 1/scissor.cut_vectors[inds,cut+1]
            values = scissor.cut_vectors[inds,cut+1]

            if sum_flows:
                sum_values.extend(values)

            if not sum_flows and len(values) != 0:
                _bins = 100
                _rng = (np.amin(values),np.max(values)*1.1)

                visualization.root_Histogram(values,
              					rng=_rng,
              					bins=_bins-int(np.sqrt(len(values))),
              					Title='{} {}'.format(flows['cut name'][cut],sample),
              					xaxis=flows['units'][cut],
              					fname='distribution_{}_{}.png'.format(cut,sample))
                
            #-----------
                
        if plot_add_info and len(scissor.survivor_inds) != 0:
            #----------- Plot additional info
            inds = np.array(scissor.survivor_inds[info_cut-1],dtype=int) # surviving indices before info_cut

            _data_x = np.array(scissor.additional_info[inds])
            _data_z = np.array(scissor.additional_info[inds])
            
            _data_x = _data_x[:,0]
            _data_x = np.delete(_data_x, np.where(_data_x == 0))
            
            if sum_flows:
                sum_data_x.extend(_data_x)            
                sum_data_z.extend(_data_z)

            _Title = 'Minimum Distance to Expected Hit of Digi Hit Before Vertex {}'.format(sample)
            _xlabel = 'Minimum Distance [cm]'
            _zlabel = 'z displacement from edge [cm]'
            _fname = 'before_dist_{}.png'.format(sample)

            det = detector.Detector()
            _xlims = [np.amin(_data_x),np.max(_data_x)*1.1]
            _zlims = [det.BoxLimits[2][0]-3000, det.BoxLimits[2][1]+3000]
            _xbins=100
            _zbins=100
            
            if not sum_flows and TwoD:
                visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, _zbins=_zbins,
    	        	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
            
            elif not sum_flows:
                visualization.root_Histogram(_data_x,
              					rng=_xlims,
              					bins=_xbins-int(np.sqrt(len(_data_x))),
              					Title=_Title,
              					xaxis=_xlabel,
              					fname=_fname,
                        logx=False,
                        logy=True)
           
           #---------------

        i += 1

    if sum_flows:
        flows['Flows fraction'] = np.round(flows['Flows'] / total_events, 8)
        
        print("Total Events: ",total_events)

        if plot_add_info:
            _xlims = [np.amin(sum_data_x),2000] #np.max(sum_data_x)*1.1]
            _zlims = [det.BoxLimits[2][0]-3000, det.BoxLimits[2][1]+3000]
            
            if TwoD:
                visualization.root_2D_Histogram(sum_data_x, sum_data_z, Title=_Title, xbins=_xbins, _zbins=_zbins,
    	        	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
                
            else:
                visualization.root_Histogram(sum_data_x,
              					rng=_xlims,
              					bins=_xbins-int(np.sqrt(len(sum_data_x))),
              					Title=_Title,
              					xaxis=_xlabel,
              					fname=_fname,
                        logx=False,
                        logy=True)
        
        if plot_cut and len(sum_values) != 0:
            _bins = 100
            _rng = (np.amin(sum_values),np.max(sum_values)*1.1)
    
            visualization.root_Histogram(sum_values,
          					rng=_rng,
          					bins=_bins-int(np.sqrt(len(values))),
          					Title='{} {}'.format(cut_name[cut],sample),
          					xaxis=cut_units[cut],
          					fname='distribution_{}_{}.png'.format(cut,sample))

    cutflow = pd.DataFrame(flows)

    print(cutflow)
    
    if len(sys.argv) == 1:
        joblib.dump(passed_events,'passed_events.joblib')


    
    
if __name__ == '__main__':

    main()



