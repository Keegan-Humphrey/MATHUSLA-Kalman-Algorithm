
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
        self.plotter = Plotter() # to be overwritten by passed plotter

    def inside_box(self, x, y, z):

        box_lims = self.det.BoxLimits
        #y_floor = self.det.LayerYLims[2][1]
        self.y_floor = self.det.y_floor

        if box_lims[0][0] < x and box_lims[0][1] > x:
           if self.y_floor < y and box_lims[1][1] > y: # only accept it if above the floor layers
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
                
            self.fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            self.fullhitindices = util.unzip(self.tree.Track_k_m_hitIndices)

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

            self.Track_floor_hits()
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
        self.cuts['Tracks'] = {'index':self.n, 'cut if':'<'} # record the index of the cut in cut vector
                                                                            # and how to operate on the value.
                                                                            # Key is the name of this function
        self.current_vector[self.n] = self.tree.NumTracks_k_m
            
            
        
    def Vertices(self):
        self.cuts['Vertices'] = {'index':self.n, 'cut if':'<'} 
        
        self.current_vector[self.n] = len(self.tree.Vertex_k_m_x)       
        

    def Fiducial_vertex(self):
            self.cuts['Fiducial_vertex'] = {'index':self.n, 'cut if':'<'} 
            
            #------------- Fiducial vertex cut
            inside = False

            for num in range(len(self.tree.Vertex_k_m_x)):
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]
                if self.inside_box(vtxx, vtxy, vtxz):
                    inside = True

            # check if at least one vertex is in the detector
            self.current_vector[3] = int(inside) # fiducial
            
            
    def Floor_hits_before_vertex(self):
            self.cuts['Floor_hits_before_vertex'] = {'index':self.n, 'cut if':'<'} 
           
            #------------- Floor hits before vertex
            # is there a hit that should be part of a track that isn't?
            
            floorveto = False

            vertex_trackIndices = self.tree.Vertex_k_m_trackIndices
            self.vertex_trackIndices = util.unzip(vertex_trackIndices)

            min_by_verts = []

            floor_wall_hits = []

            for hit in range(len(self.tree.Digi_y)): # loop over the hits
                if self.in_layer(self.tree.Digi_y[hit]) <= 2: 
                    floor_wall_hits.append(hit)
                
                elif self.tree.Digi_z[hit] <= self.det.z_wall: # check if layer index is a floor layer index or in the wall
                    floor_wall_hits.append(hit)
                    

            for num in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vertexveto = False
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]

                if not self.inside_box(vtxx, vtxy, vtxz): # make sure the vertex is fiducial (above the floor)
                    vertexveto = True # veto the vertex if its outside the detector 

                
                if not vertexveto:
                    ''' x-z distance of track '''
                    min_dist = 1e6 # arbitrary, large value to start us off
              
                    floor_wall_hits_before_vertex = []
              
                    for hit in floor_wall_hits:
                        if self.tree.Vertex_k_m_t[num] + 2 * self.tree.Vertex_k_m_ErrorT[num] > self.tree.Digi_time[hit]: # check if the floor or wall hit happened before the vertex
                            floor_wall_hits_before_vertex.append(hit)
                            
                    for ind in self.vertex_trackIndices[num]: # Determine where the tracks in the vertex are expected to hit the wall and floor
                        ind = int(ind)
    
                        x = self.tree.Track_k_m_x0[ind]
                        y = self.tree.Track_k_m_y0[ind]
                        z = self.tree.Track_k_m_z0[ind]
    
                        vx = self.tree.Track_k_m_velX[ind]
                        vy = self.tree.Track_k_m_velY[ind]
                        vz = self.tree.Track_k_m_velZ[ind]
                        
                        try:
                            del_t = (self.y_floor - y) / vy
                            
                            #expected_x_floor.append(x + vx * del_t)
                            #expected_z_floor.append(z + vz * del_t)
                            
                            expected_x_floor = x + vx * del_t
                            expected_z_floor = z + vz * del_t
                            
                            del_t = (self.det.z_wall - z) / vz
                            
                            expected_x_wall = x + vx * del_t
                            expected_y_wall = y + vy * del_t
                        
                        except ZeroDivisionError:
                            continue

                        for hit in floor_wall_hits_before_vertex:
                            min_dist_floor = np.amin(np.sqrt((expected_x_floor - self.tree.Digi_x[hit])**2 + (expected_z_floor - self.tree.Digi_z[hit])**2)) # min dist to a track for the hit in plane of the floor
                            min_dist_wall = np.amin(np.sqrt((expected_x_wall - self.tree.Digi_x[hit])**2 + (expected_y_wall - self.tree.Digi_y[hit])**2)) # min dist to a track for the hit in plane of the wall
                            
                            min_dist = np.amin([min_dist, min_dist_floor, min_dist_wall])


                    if min_dist < 1500: # [cm] Make sure closest hit to the track is within 3 m to veto
                                                 # ******** need to find a slick way to let this value be cut on independently

		                            # *** look at hit time as well as distance -> pass to plotter
                            vertexveto = True
                        

                if not vertexveto: # There's a vertex in the event that isn't vetoed
                    break

            if len(self.tree.Vertex_k_m_x) != 0:
                if vertexveto:
                    floorveto = True

            self.current_vector[self.n] = int(not floorveto) # cut if floorveto == True => 0 == current_vector[4] (< 1 is true)

                

    def Track_floor_hits(self):
            self.cuts['Track_floor_hits'] = {'index':self.n, 'cut if':'True'}
            
            #---------- Reject tracks with digi hits in the floor
            
            track_floor_veto = True # if we make it through the loop unscathed, veto the event
            
            for vert in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vertex_veto = False    
                
                for trk in self.vertex_trackIndices[vert]: # loop over tracks in the vertex 
                    trk = int(trk)
                    if self.tree.Track_k_m_y0[trk] <= self.y_floor: # check if the lowest hit is in the floor
                        vertex_veto = True
                        break
                        
                        #***** add wall veto as well
                
                if not vertex_veto: # if there is a vertex with without floor hits, don't veto the event
                    track_floor_veto = False
                    break
            
            self.current_vector[self.n] = track_floor_veto


    def Opening_angle(self):
            self.cuts['Opening_angle'] = {'index':self.n, 'cut if':'<'}
    
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
            self.cuts['Topological'] = {'index':self.n, 'cut if':'<'}

            #---------------- check if vertices are below track hits within some sigma (Topological Cut)
            topoveto = 0

            #self.fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            #self.fullhitindices = util.unzip(self.tree.Track_k_m_hitIndices)

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
            self.cuts['Vertex_track_beta'] = {'index':self.n, 'cut if':'<'}
             
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
            self.cuts['Track_hit_diffs'] = {'index':self.n, 'cut if':'<'}
             
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
            self.cuts['Exp_hits'] = {'index':self.n, 'cut if':'<'}
             
            #------------- Expected Hits Cut
            # if there are no hits in the floor / wall, could it be that the particle avoided them?

            centre_floor = [(self.det.BoxLimits[0][0] + self.det.BoxLimits[0][1]) / 2, (self.det.BoxLimits[2][0] + self.det.BoxLimits[2][1]) / 2] # [x,z] [cm,cm]
            centre_wall = [centre_floor[0], 6000.0 + 2000.0 / 2.0] # [x,y] [cm,cm] see globals.hh for explanation on y

            floor_ys = [(self.det.LayerYLims[i][0] + self.det.LayerYLims[i][1]) / 2 for i in range(3)]
            floor_y = floor_ys[2]
  
            '''
            hits_in_floor = False # see if there are any hits in the floor *** need to change this to hits nearby the track

            for i in range(3):
                if np.any(np.array(self.tree.Digi_y) == np.sum(self.det.LayerYLims[i]) / 2):
                    hits_in_floor = True
                    
            hits_in_wall = False # see if there are any hits in the wall
            
            if np.any(np.array(self.tree.Digi_z) <= self.det.z_wall):
                hits_in_wall = True
            '''       
            hits_in_wall = False # see if there are any hits in the wall
            hits_in_floor = False # see if there are any hits in the floor *** need to change this to hits nearby the track

            for hit in range(len(self.tree.Digi_y)): # loop over the hits
                if self.in_layer(self.tree.Digi_y[hit]) <= 2:
                    hits_in_floor = True
                    break
                
                if self.tree.Digi_z[hit] <= self.det.z_wall:
                    hits_in_wall = True
                    break
                
            max_centre_dist = -1e6
            second_max_centre_dist = -1e6 # large and negative arbitrary value to start us off
                    
            if len(self.vertex_trackIndices) > 0 and not (hits_in_floor or hits_in_wall): # if we have vertices made but don't find any 
                                                                                          # hits that would let us reject the event in 
                                                                                          # the hits before vertex cut, and it's background,
                                                                                          # we expect incorrect track reconstruction to be the cause
                                                                                          # ie. tracks will be close to the edge of wall / floor
                for num in range(len(self.tree.Vertex_k_m_x)):
                    for ind in self.vertex_trackIndices[num]: # Determine where the tracks in the vertex are expected to hit the wall and floor
                        ind = int(ind)
        
                        x = self.tree.Track_k_m_x0[ind]
                        y = self.tree.Track_k_m_y0[ind]
                        z = self.tree.Track_k_m_z0[ind]
        
                        vx = self.tree.Track_k_m_velX[ind]
                        vy = self.tree.Track_k_m_velY[ind]
                        vz = self.tree.Track_k_m_velZ[ind]
                        
                        try:
                            del_t = (self.y_floor - y) / vy
                            
                            expected_x_floor = x + vx * del_t
                            expected_z_floor = z + vz * del_t
                            
                            del_t = (self.det.z_wall - z) / vz
                            
                            expected_x_wall = x + vx * del_t
                            expected_y_wall = y + vy * del_t
                        
                        except ZeroDivisionError:
                            continue
                        
                        edge_dist_x_floor = (self.det.BoxLimits[0][1] - self.det.BoxLimits[0][0]) / 2 - np.abs(expected_x_floor - centre_floor[0]) # distance from expected hits
                        edge_dist_z_floor = (self.det.BoxLimits[2][1] - self.det.BoxLimits[2][0]) / 2 - np.abs(expected_z_floor - centre_floor[1]) # to edge of floor for the vertex
                        edge_dist_x_wall = (self.det.WallLimits[0][1] - self.det.WallLimits[0][0]) / 2 - np.abs(expected_x_wall - centre_wall[0]) # distance from expected hits
                        edge_dist_y_wall = (self.det.WallLimits[1][1] - self.det.WallLimits[1][0]) / 2 - np.abs(expected_y_wall - centre_wall[1]) # to edge of wall for the vertex
                        
                        max_dist = np.max([edge_dist_x_wall, edge_dist_y_wall,edge_dist_x_floor, edge_dist_z_floor]) # take the second closest track to the centre of sensitive layer in any direction
                                                                                                                     # if it passes cut, then there are two tracks in the vertex that pass the cut
                                                
                        if max_dist > max_centre_dist:
                            second_max_centre_dist = max_centre_dist
                            max_centre_dist = max_dist
                            
                        elif max_dist > second_max_centre_dist:
                            second_max_centre_dist = max_dist
                        
                
                self.current_vector[self.n] = max_centre_dist                                         
                
                 
            if (hits_in_floor or hits_in_wall):
                self.current_vector[self.n] = 5000.0 # [cm] if there's hits in the wall or floor, assign largest possible value (at detector centre) so it doesn't get cut
            
                  

    def No_floor_hits(self):
        self.cuts['No_floor_hits'] = {'index':self.n, 'cut if':'<'}
             
        ''' Veto the event if it has floor hits '''
    
        floor_veto = False
    
        for hit in range(len(self.tree.Digi_y)): # loop over the hits
            if self.in_layer(self.tree.Digi_y[hit]) <= 2: # check if layer index is a floor layer index
                floor_veto = True
                break

        if floor_veto:
            self.current_vector[self.n] = True  
            

    def Missing_hit_sum(self):
        self.cuts['Missing_hit_sum'] = {'index':self.n, 'cut if':'<'}
             
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
        self.cuts['Chi_ndof_cut'] = {'index':self.n, 'cut if':'<'}
    
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
        self.plotter = Plotter() # to be overwritten by passed plotter


    def store_space(self):

        space = sample_space(self.file)
        space.plotter = self.plotter
        
        self.cut_vectors = space.get_states_m()
        self.func_dicts = space.cuts
        
        self.plotter = space.plotter # update the plotter to collect information from the parameter space
      
        self.additional_info = space.additional_info
      
        self.events = len(self.cut_vectors)
        
        print('number of events is: ',len(self.cut_vectors))


    def cut_dict(self, cut_options, permutation):

        global ncuts

        self.flows = np.zeros(ncuts)
        local_vectors = np.copy(self.cut_vectors)
        self.survivor_inds = []
        
        if local_vectors.size != 0:
            
            for i in permutation: # do the cuts in the order specified by permutation
            
                opt_dict = cut_options[str(i)]
                func = opt_dict['func_name']
                func_dict = self.func_dicts[func]
                
                n = func_dict['index'] # cut vector index for this cut
                
                if opt_dict['on?'] == 1: 
                    cut_parameter = opt_dict['cut parameter']
                
                    if func_dict['cut if'] == '<':
                        inds = np.where(local_vectors[:,n] < cut_parameter) 
                        
                    elif func_dict['cut if'] == 'True':
                        inds = np.where(local_vectors[:,n] == True)   
                        
                    # add an elif to create a new cut condition
                        
                else:
                    inds = [] # don't cut any events
                     
                     
                local_vectors = np.delete(local_vectors, inds, 0)

                self.flows[i] += np.shape(local_vectors)[0]
                
                self.survivor_inds.append(local_vectors[:,0])

        self.survivors = self.flows[-1]
        
        
        
class Plotter:

    def __init__(self):
        ''' Assign to attributes to collect data for plots '''
        
        self.min_dist_floor = []
        self.min_dist_wall = []
        self.sample = 'sim data' # type of sample being used to make plots
       
        
    def merge_a_plotter(self, pltr):
        
        self.min_dist_floor.extend(pltr.min_dist_floor)
        self.min_dist_wall.extend(pltr.min_dist_wall) # change data handling to dictionary and loop over keys here
       
        
    def plot(self):
        ''' Generate the plots you want with any of the methods! '''
    
        self.plot_min_exp_dist()
    
        
    def plot_min_exp_dist(self):
    
        _Title = 'Minimum Distances to Expected Hit of Digi Hit Before Vertex {}'.format(self.sample)
        _xlabel = 'Minimum Distance to a floor hit [cm]'
        _zlabel = 'minimum distance to a wall hit [cm]'
        _fname = 'min_dist_{}.png'.format(self.sample)

        _data_x = self.min_dist_floor
        _data_z = self.min_dist_wall

        det = detector.Detector()
        _xlims = [np.amin(_data_x),np.max(_data_x)*1.1]
        _zlims = [np.amin(_data_z),np.max(_data_z)*1.1]
        _xbins=100
        _zbins=100
        
        visualization.root_Histogram(_data_z, rng=None, ft_rng=None, bins=0, Title=_Title+" (wall)", xaxis=_zlabel, logx=False, logy=False, fname='min_dist_wall_{}.png'.format(self.sample))
        visualization.root_Histogram(_data_x, rng=None, ft_rng=None, bins=0, Title=_Title+" (floor)", xaxis=_xlabel, logx=False, logy=False, fname='min_dist_floor_{}.png'.format(self.sample))
        
        visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, zbins=_zbins,
	      	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
        
        

def main():
  
    if len(sys.argv) > 1:
        # default set for job submission with files and write locations passed as arguments
        
        if len(sys.argv) != 4:
            print('Usage: python cutflow.py <tracker tree> <save directory> <integer identifier>')
            return
    
        plot_cut = False
        plot_add_info = False
        plot_obj = False
        
        sum_flows = True # True <=> Background
        load = False
        save = True
        TwoD = False
        
        files = [sys.argv[1]]

    else:
        # plotting booleans (to be deprecated soon!)
        plot_cut = False
        plot_add_info = False
        plot_obj = False
        
        sum_flows = True # True <=> Background / sum over data in files for flows ***** need to adress sum_flows or load booleans in below code
        load = False
        save = False
        TwoD = False # whether additional info is 2D
        
        if load:
            #load_files_dir = "/home/keeganh/scratch/job_test/W_sample_dir/run6/analysis_data/21_01_22/"
            
            load_files_dir = "/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run3/analysis_data/30_01_22/"
            files = [filename for filename in glob.iglob(load_files_dir+'/**/scissor_*.joblib', recursive=True)]
            
            #load_files_dir = "/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/analysis/save_files/"
            #files = [load_files_dir+'scissor_{}.joblib'.format(sample) for sample in ['h10_0','h2_1','qq_2']]
            
        else:
            if sum_flows: # Background
                directory_3 = '/home/keeganh/scratch/job_test/W_sample_dir/run3/18_12_21/11_24_04/trees/'
                #directory_4 = '/home/keeganh/scratch/job_test/W_sample_dir/run4/09_01_22/09_11_57/trees/'
                
                files = [filename for filename in glob.iglob(directory_3+'stat_*.root', recursive=True)]
                #files.extend([filename for filename in glob.iglob(directory_4+'stat_*.root', recursive=True)])
        
            else: # signal
                files = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_0_0.root',
                    '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_2_0.root']
                    #'/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm//08_01_22/17_50_20/trees/stat_0_0.root']
        #            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_4_0.root']
           
    if len(files) == 0:
        print("I need at least 1 file to run!")
        return

    cut = 9 # which cut to plot (index in cut_options)
    info_cut = 3 # which cut we take additional info from 
    
    sum_values = []
    sum_data_x = []
    sum_data_z = []
    drawers = [] # to hold the scissors

    flows = dict()
    passed_events = dict()
    
    option = ['cut name', 'cut parameter', 'units', 'on?'] # [name, value of the cut, units of what is being cut on, whether the cut is used (1 => yes)]
    func_name = 'func_name'

    # specify cuts and parameters
    cut_options = {'0' :{option[0]:'2 Tracks'                     ,option[1]:2       ,option[2]:'tracks'      ,option[3]:1 , func_name:'Tracks' },
                   '1' :{option[0]:'Vertices'                     ,option[1]:1       ,option[2]:'verts'       ,option[3]:1 , func_name:'Vertices' },
                   '2' :{option[0]:'Fiducial Vertex'              ,option[1]:1       ,option[2]:'bool'        ,option[3]:1 , func_name:'Fiducial_vertex' },
                   '3' :{option[0]:'Floor/Wall Hits Before Vertex',option[1]:1       ,option[2]:'bool'        ,option[3]:1 , func_name:'Floor_hits_before_vertex' },
                   '4' :{option[0]:'No track hits in the floor'   ,option[1]:True    ,option[2]:'bool'        ,option[3]:1 , func_name:'Track_floor_hits' },
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
                                                    
    cut = np.where(np.array(permutation) == cut)[0][0] # *** find a slicker way to handle cut plotting -> encorporate into plotter?
    
    for key in flows.keys():
        flows[key] = flows[key].take(permutation,0) # permute each entry in flows dictionary


    if sum_flows and plot_obj:
        pltr = Plotter() # object to collect and organize information in the cuts to be plotted
 
    total_events = 0

    i = 0

    for file in files:

        if sum_flows or len(sys.argv) > 1:
            sample = 'W'

        else:
            sample = ['h10','h2','qq'][i] # particle type for each sample in files 


        if not load:
            print(file)
        
            scissor = scissors(file)
            
            if plot_obj:
                if not sum_flows:
                    pltr = Plotter()
                
                pltr.sample = sample
                
                scissor.plotter = pltr

            try:
                scissor.store_space() # read the trees and gather data for cutting from each event
            
            except ReferenceError: # tree is empty
                print("Sorry, that tree doesn't have anything in it")
                i += 1
                continue

            pltr = scissor.plotter # update the global plotter opject to include info from the current parameter space
            
            
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

        if plot_obj and not sum_flows:
            pltr.plot() # generate plots for this sample with plotter
            
        if load:
            if i % 100 == 0:
                print(i," files processed")
            
            try:
                scissor = joblib.load(file)
                #scissor = joblib.load('save_files/scissor_{}_{}.joblib'.format(sample,i)) # load the scissor objects so we don't
                                                                                          # need to read and process the data again
                                                                                          # (just need to cut with chosen parameters)                                                              
                if plot_obj:
                    pltr.sample = sample
                    
                    pltr.merge_a_plotter(scissor.plotter)
                    
                
            except (FileNotFoundError, EOFError, OSError) as error:
                print("Sorry, I encountered an error with that file")
                #print(error)
                i += 1
                continue
                
        scissor.cut_dict(cut_options, permutation)
        
        #passed_events[file] = scissor.survivor_inds
        passed_events[scissor.file] = scissor.survivor_inds
        
        
        file_flows = scissor.flows.astype(int)
        file_flows = file_flows.take(permutation,0) # permute the flows to the order the cuts were performed
        
        if sum_flows:
            try:
                flows['Flows'] += file_flows
                             
            except:
                flows['Flows'] = file_flows

            total_events += scissor.events

        else:
            flows[sample] = file_flows
            flows[sample+' fraction'] = np.round(file_flows / scissor.events, 4)
            

        if plot_cut and len(scissor.survivor_inds) != 0:
            #------- Plot the Chosen cut
            inds = np.array(scissor.survivor_inds[cut-1],dtype=int) # surviving indices before cut

#            values = 1/scissor.cut_vectors[inds,cut+1]
            #values = scissor.cut_vectors[inds,cut+1]
            values = scissor.cut_vectors[inds,scissor.func_dicts[cut_options[str(permutation[cut])]['func_name']]['index']]

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

        if plot_obj:
            pltr.plot() # generate plots for the run with plotter


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
          					Title='{} {}'.format(flows['cut name'][cut],sample),
          					xaxis=flows['units'][cut],
          					fname='distribution_{}_{}.png'.format(cut,sample),
                    logy=True)


    cutflow = pd.DataFrame(flows)

    print(cutflow)
    
    joblib.dump(cutflow,'flows.joblib')
    
    if len(sys.argv) == 1 and save:
        joblib.dump(passed_events,'passed_events.joblib')


    
    
if __name__ == '__main__':

    main()



