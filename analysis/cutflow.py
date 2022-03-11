
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

ncuts = 18

events_passed = 0
ev_w_veto_hits = 0

class sample_space():


    def __init__(self,file):

        self.file = file
        
        self.cuts = dict() # information about how to use the cuts
        self.plotter = Plotter()
        
        self.det = detector.Detector()
        self.y_floor = self.det.y_floor
        self.z_wall = self.det.z_wall
        
        self.IP_to_wall_slope = (self.z_wall / self.det.WallLimits[1][1]) # z / y to top of wall from IP
        
       
    def inside_box(self, x, y, z):

        box_lims = self.det.BoxLimits
        
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

        self.tree.SetBranchStatus("Track_k_m_velX", 1)
        self.tree.SetBranchStatus("Track_k_m_velY", 1)
        self.tree.SetBranchStatus("Track_k_m_velZ", 1)
        self.tree.SetBranchStatus("Track_k_m_hitIndices", 1)
        self.tree.SetBranchStatus("Track_k_m_expected_hit_layer", 1)

        self.tree.SetBranchStatus("Track_k_m_t0", 1)
        self.tree.SetBranchStatus("Track_k_m_x0", 1)
        self.tree.SetBranchStatus("Track_k_m_y0", 1)
        self.tree.SetBranchStatus("Track_k_m_z0", 1)
        self.tree.SetBranchStatus("Track_k_smooth_chi_sum",1)

        cut_vectors = []
        
        self.additional_info = np.zeros((int(self.tree.GetEntries()),2)) # Fill me up with info if you want internal 
                                                                         # access to the cuts to plot something new
        for event_number in range(int(self.tree.GetEntries())):

            self.event_number = event_number
            
            self.event_info = event_info(event_number)

            self.current_vector = np.zeros(1+ncuts) # first will be the event number
            self.n = 0

            self.tree.GetEntry(event_number)

            if event_number % 1000 == 0:
                print("event:", event_number)
                pass
                
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
            
            #------------- cut if the vertex is outside the box
            
            self.Fiducial_vertex()
            self.n += 1
            
            #------------- Optional Cuts
            
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
            
            self.Hits_per_track()
            self.n += 1
            
            self.Fiducial_leniency()
            self.n += 1

            self.Delta_ray_cut()
            self.n += 1
            
            self.Fiducial_vertex_time()
            self.n += 1

            self.Chi_ndof_cut()
            self.n += 1
            
            #------------ All done
            cut_vectors.append(self.current_vector)
            
            self.plotter.event_infos.append(self.event_info)

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

            vtxx, vtxy, vtxz = 0, 0, 0

            for num in range(len(self.tree.Vertex_k_m_x)):
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]
                if self.inside_box(vtxx, vtxy, vtxz):
                    inside = True
                    
                    #self.event_info.data_dict['vert pos'].append([vtxx,vtxy,vtxz]) # add all vertices
                    
            # check if at least one vertex is in the detector
            self.current_vector[3] = int(inside) # fiducial
            
            #self.plotter.data_dict['vert pos'].append([vtxx,vtxy,vtxz]) # add one representative vertex or 0s if none

     
    def Floor_hits_before_vertex(self):
            self.cuts['Floor_hits_before_vertex'] = {'index':self.n, 'cut if':'<'} 
           
            #------------- Floor hits before vertex
            # Is there a hit that should be part of a track that isn't?
            #
            # keep track of the smallest distance from the expected hit position of a track in a vertex
            # to a hit in the floor or wall, if this is smaller than the cut value, we think the particle that 
            # formed the vertex could have been a W and the tracker simply missed a floor / wall hit
            #
            # Essentially we add them back in by hand and veto the event because it has a track with floor hits
            
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
                    
            min_dist = 1e6 # arbitrary, large value to start us off
            min_time = 1e6 # also keep track of time for plots

            for num in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vertexveto = False
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]

                if not self.inside_box(vtxx, vtxy, vtxz): # make sure the vertex is fiducial (covered by wall / floor and in the detector)
                    continue
                
                if not vertexveto:
                    ''' x-z distance of track '''
                    
                    floor_wall_hits_before_vertex = []
              
                    for hit in floor_wall_hits:
                        if self.tree.Vertex_k_m_t[num] + 2 * self.tree.Vertex_k_m_ErrorT[num] > self.tree.Digi_time[hit]: # check if the floor or wall hit happened before the vertex
                            floor_wall_hits_before_vertex.append(hit)
                            
                    for ind in self.vertex_trackIndices[num]: # Determine where the tracks in the vertex are expected to hit the wall and floor
                        ind = int(ind)
    
                        x = self.tree.Track_k_m_x0[ind]
                        y = self.tree.Track_k_m_y0[ind]
                        z = self.tree.Track_k_m_z0[ind]
                        t = self.tree.Track_k_m_t0[ind]
    
                        vx = self.tree.Track_k_m_velX[ind]
                        vy = self.tree.Track_k_m_velY[ind]
                        vz = self.tree.Track_k_m_velZ[ind]
                        
                        try:
                            del_t_floor = (self.y_floor - y) / vy
                            
                            expected_x_floor = x + vx * del_t_floor
                            expected_z_floor = z + vz * del_t_floor
                            
                            del_t_wall = (self.det.z_wall - z) / vz
                            
                            expected_x_wall = x + vx * del_t_wall
                            expected_y_wall = y + vy * del_t_wall
                        
                        except ZeroDivisionError:
                            continue

                        dist_floors, dist_walls = [], []
                        dt_floors, dt_walls = [], []
                        
                        for hit in floor_wall_hits_before_vertex:
                            dist_floor = np.sqrt((expected_x_floor - self.tree.Digi_x[hit])**2 + (expected_z_floor - self.tree.Digi_z[hit])**2) # min dist to a track for the hit in plane of the floor
                            dist_wall = np.sqrt((expected_x_wall - self.tree.Digi_x[hit])**2 + (expected_y_wall - self.tree.Digi_y[hit])**2) # min dist to a track for the hit in plane of the wall
                            
                            dt_floor = abs(self.tree.Digi_time[hit] - (t + del_t_floor))
                            dt_wall = abs(self.tree.Digi_time[hit] - (t + del_t_wall))
                            
                            min_dist = np.amin([min_dist, dist_floor, dist_wall])
                            min_time = np.amin([min_time, dt_floor, dt_wall])
                    
            self.current_vector[self.n] = min_dist 
            
            ## add timing info to plotter ******
                

    def Track_floor_hits(self):
            self.cuts['Track_floor_hits'] = {'index':self.n, 'cut if':'True'}
    
            #--------- Reject events with digi hits in the floor (by event)
            
            track_floor_veto = False
            
            for trk in range(len(self.tree.Track_k_m_y0)): # loop over tracks in the event 
                trk = int(trk)
                
                for hit in self.fullhitindices[trk]:
                    if self.tree.Digi_y[hit] <= self.y_floor: # check if the hit is in the floor
                        track_floor_veto = True
                        break
                        
                    elif self.tree.Digi_z[hit] <= self.z_wall: # check if the hit is in the wall
                        track_floor_veto = True
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
                                                                                          # the floor_wall_hits_before_vertx, and it's a W background,
                                                                                          # we expect incorrect track reconstruction velocity to be the cause
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
                            
                            #expected_x_floor = x + vx * del_t
                            expected_z_floor = z + vz * del_t
                            
                            del_t = (self.det.z_wall - z) / vz
                            
                            #expected_x_wall = x + vx * del_t
                            expected_y_wall = y + vy * del_t
                        
                        except ZeroDivisionError:
                            continue
                        
                        # calculate the distance from the edge of sensitive layers (floor or wall)
                        # symmetries modded out so < 0 => outside the detector in each direction
                        
                        if expected_y_wall >= self.y_floor: # doesn't go through the floor, we prefer this one!
                            #edge_dist_y_wall = (self.det.WallLimits[1][1] - self.det.WallLimits[1][0]) / 2 - np.abs(expected_y_wall - centre_wall[1]) # to edge of wall for the vertex
                            edge_dist_y_wall = (self.det.WallLimits[1][1] - self.det.WallLimits[1][0]) / 2 - (centre_wall[1] - expected_y_wall) # to edge of wall for the vertex
                            #edge_dist_x_wall = (self.det.WallLimits[0][1] - self.det.WallLimits[0][0]) / 2 - np.abs(expected_x_wall - centre_wall[0]) # distance from expected hits
                        
                            #self.plotter.data_dict['Exp Pos x wall'].append(edge_dist_x_wall)
                            #self.plotter.data_dict['Exp Pos y wall'].append(edge_dist_y_wall)
                        
                            #min_dist = np.amin([edge_dist_x_wall, edge_dist_y_wall])
                            min_dist = edge_dist_y_wall # remove x and positive edge cuts -> just corner between wall and floor
                        
                        else: # doesn't doesn't go through the wall, we prefer this one!
                            #edge_dist_z_floor = (self.det.BoxLimits[2][1] - self.det.BoxLimits[2][0]) / 2 - np.abs(expected_z_floor - centre_floor[1]) # to edge of floor for the vertex
                            edge_dist_z_floor = (self.det.BoxLimits[2][1] - self.det.BoxLimits[2][0]) / 2 - (centre_floor[1] - expected_z_floor) # to edge of floor for the vertex
                            #edge_dist_x_floor = (self.det.BoxLimits[0][1] - self.det.BoxLimits[0][0]) / 2 - np.abs(expected_x_floor - centre_floor[0]) # distance from expected hits
                        
                            #self.plotter.data_dict['Exp Pos x floor'].append(edge_dist_x_floor)
                            #self.plotter.data_dict['Exp Pos z floor'].append(edge_dist_z_floor)
                        
                            #min_dist = np.amin([edge_dist_z_floor, edge_dist_x_floor]) # **** these are from either edge, we should really have distance from corner of wall and floor
                                                                                        # ie. only look at + z and + y directions
                            min_dist = edge_dist_z_floor # remove x and positive edge cuts -> just corner between wall and floor
                          
                        # if the distance from the edge is larger (closer to det centre) than other tracks found,
                        # keep track of it 
                        if min_dist > max_centre_dist:
                            second_max_centre_dist = max_centre_dist
                            max_centre_dist = min_dist
                            
                        # keep track of second closest to the centre; we will cut on this to make sure 2 tracks pass
                        elif min_dist > second_max_centre_dist:
                            second_max_centre_dist = min_dist
                
                self.current_vector[self.n] = second_max_centre_dist                                        
                
                 
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
        
        # *** add condition on number of tracks / vertices to increase quark efficiency
        
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
        

    def Hits_per_track(self):
        ''' make sure at least two tracks in a vertex have >= ___ # of hits '''
        
        self.cuts['Hits_per_track'] = {'index':self.n, 'cut if':'<'}
        
        track_hits = []
        
        for vert in range(len(self.vertex_trackIndices)): # vertex
            vert = int(vert)
            
            vert_tracks = []
            
            for trk in self.vertex_trackIndices[vert]: # track
                trk = int(trk)
                
                vert_tracks.append(abs(len(self.fullhitindices[trk])))
                
            vert_tracks = np.abs(vert_tracks)
            #vert_tracks.sort(reverse=True)
            vert_tracks = np.sort(vert_tracks)
            
            track_hits.append(vert_tracks[-2]) # add the 2nd highest number of hits per track in the vertex
                                                # so two tracks are required to make the cut
            
        if len(track_hits) != 0:
            self.current_vector[self.n] = np.max(track_hits) # choose the lowest chis
            
    
    def Fiducial_leniency(self):
    
        self.cuts['Fiducial_leniency'] = {'index':self.n, 'cut if':'<'}
    
        # Could the particle have come in above the wall?
        #
        # if the vertex is fiducial, calculate the distance to the line through the top of the wall and IP 
        #
        # distance^2 to a line in the plane through the origin (the IP) with slope m from a point (z, y) is 
        # ( z - y m ) / ( 1 + m^2 ) , with m = self.IP_to_wall_slope for us
    
        vtxx, vtxy, vtxz = 0, 0, 0
        
        max_dist_to_slope = -1e5

        for num in range(len(self.tree.Vertex_k_m_x)):
            vtxx, vtxy, vtxz,  = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]
             
            if not self.inside_box(vtxx, vtxy, vtxz):
                continue
                
            else:
                vtx_dist_to_slope = ( vtxz - vtxy * self.IP_to_wall_slope ) / ( 1 + self.IP_to_wall_slope**2 )
                max_dist_to_slope = max(max_dist_to_slope, vtx_dist_to_slope)
        
        self.event_info.data_dict['slope dist'].append(max_dist_to_slope)
        
        self.current_vector[self.n] = max_dist_to_slope
        
        
    def Delta_ray_cut(self):
        
        self.cuts['Delta_ray_cut'] = {'index':self.n, 'cut if':'<'}
        
        '''
        do a rough pairwise tracking run 
        with a hit in the floor / wall before tracking hits and unused hits in tracking layers. We then keep 
        pairs that would have 0.8 < beta < 1.2 (like the tracker) and use that to compute the closest approach distance with made tracks. 
        '''
        closest_approach = 1e6 # start us off large and arbitrary
        lowest_pair = None
        
        tell_me_baby_whats_your_story = False
        
        vertices = self.current_vector[self.cuts['Vertices']['index']]
        fiducial_vertex = bool(self.current_vector[self.cuts['Fiducial_vertex']['index']])
        #no_fiducial_vertex = False
        no_track_hits_in_floor = bool(self.current_vector[self.cuts['Track_floor_hits']['index']])
        
        if vertices > 0 and fiducial_vertex and not no_track_hits_in_floor: # don't compute anything if there are no vertices
            #self.current_vector[self.n] = closest_approach
            #return
            
            '''
            global events_passed                
            events_passed += 1
            print(events_passed,' Events passed vetos ---------')
            '''
        
            # find all hits not involved in a track
            track_hit_indices = set([hit for track in self.fullhitindices for hit in track])
            all_hit_indices = set(np.arange(len(self.tree.Digi_y)))
            
            not_track_hits = all_hit_indices - track_hit_indices
            
            # divide hits into tracking vs floor / wall hits 
            veto_hits = [] # floor / wall hits (not in a track)
            tracking_hits = [] # tracking layer hits (not in a track)
    
            for hit in not_track_hits: # loop over the hits not in a track
                hit = int(hit)
                
                if self.tree.Digi_z[hit] <= self.det.z_wall:
                    veto_hits.append(hit) # its in the wall
        
                else:
                    if self.in_layer(self.tree.Digi_y[hit]) <= 2:
                        veto_hits.append(hit) # its in the floor
                        
                    else:
                        tracking_hits.append(hit) # its above the floor => tracking layers
        
            
            if len(veto_hits) != 0:
                # find time of earliest tracking layer hit not in a track

                '''
                global ev_w_veto_hits 
                ev_w_veto_hits += 1
                print(ev_w_veto_hits, ' have veto hits')
                print('event ',self.event_number,' I have {} veto hits'.format(len(veto_hits)))
                '''
                
                min_vertex_t = 1e6 # [ns] time of earliest tracking layer hits
                min_tracking_t = 1e6 # [ns] time of earliest tracking layer hits
                
                for vertex in range(len(self.tree.Vertex_k_m_t)): # loop over lists of vertices
                    min_vertex_t = min(min_vertex_t, self.tree.Vertex_k_m_t[vertex])
                ''''''
                
                # **** switch back to the version below (tracking layer hits instead of vertices for min time)
                '''    '''    
                for hit in tracking_hits:
                    min_tracking_t = min(min_tracking_t, self.tree.Digi_time[hit])
                
                if tell_me_baby_whats_your_story:
                    print('min vertex t',min_vertex_t)
                    print('min tracking t ',min_tracking_t)
                
                min_vertex_t = max(min_vertex_t, min_tracking_t)
                
                # find all the floor / wall hits before all the tracking layer hits (ie. not backscatters)
                early_veto_hits = []
               
                if tell_me_baby_whats_your_story:
                    print('looking for veto hits')
                
                for hit in veto_hits:
                    if self.tree.Digi_time[hit] < min_vertex_t:
                        
                        if tell_me_baby_whats_your_story:
                            print(self.tree.Digi_time[hit],' early veto hit')
                        
                        early_veto_hits.append(hit)
               
                if tell_me_baby_whats_your_story:
                    print(len(early_veto_hits),' early veto hits')
                
                # find all pairwise combinations of veto / tracking hits that have beta ~ c (0.8 < beta < 1.2 like real tracks)
                paired_tracks = []
                
                for v_hit in early_veto_hits: 
                    for t_hit in tracking_hits:
                    
                        pair = pair_track(v_hit, t_hit)
                        
                        pair.get_st_points(self.tree)
                        pair.find_beta()
                        
                        if tell_me_baby_whats_your_story:
                            print(t_hit,' beta is ',pair.beta)
                        
                        if 0.8 < pair.beta < 1.2:
                            paired_tracks.append(pair)
                
                # find smallest distance between paired_tracks and tracks in a vertex at any time (closest approach distance, see docs)  
                
                for pair in paired_tracks:
                    # no need anymore, done above already (called in find_beta())
                    #pair.find_velocity()
                
                    for vertex in self.fulltrackindices: # loop over lists of vertices
                        for track in vertex: # loop over tracks in each vertex list
                            track = int(track)
                            
                            pair.find_closest_approach(self.tree, track)
                            
                            if tell_me_baby_whats_your_story:
                                print('Closest approach ',pair.closest_approach)
                            
                            #closest_approach = min(closest_approach, pair.closest_approach)
                            if pair.closest_approach < closest_approach:
                                closest_approach = pair.closest_approach
                                lowest_pair = pair
                
        if lowest_pair != None:
            self.event_info.data_dict['pair reco beta'].append(lowest_pair.beta)
            self.event_info.data_dict['pair closest approach'].append(lowest_pair.closest_approach)
            
            if tell_me_baby_whats_your_story:
                print('lowest closest approach ',lowest_pair.closest_approach)
            
        else:
            self.event_info.data_dict['pair reco beta'].append(1e6)
            self.event_info.data_dict['pair closest approach'].append(1e6)
            
        self.current_vector[self.n] = closest_approach
        

    def Fiducial_vertex_time(self):
        self.cuts['Fiducial_vertex_time'] = {'index':self.n, 'cut if':'True'}
    
        # Remove mis timed vertices, if all vertices happen before all the floor hits
        # veto the event, it can avoid our cuts and there is no valid physical way signal
        # can manage this
        
        all_hit_indices = np.arange(len(self.tree.Digi_y))
        
        max_vertex_t = -1e6 # [ns] time of latest vertex hit
        
        for vertex in range(len(self.tree.Vertex_k_m_t)): # loop over vertices
            max_vertex_t = max(max_vertex_t, self.tree.Vertex_k_m_t[vertex])
        
        veto_hits = [] # floor / wall hits
            
        for hit in all_hit_indices:
            hit = int(hit)
            
            if self.tree.Digi_z[hit] <= self.det.z_wall:
                veto_hits.append(hit) # its in the wall
                
            elif self.in_layer(self.tree.Digi_y[hit]) <= 2:
                veto_hits.append(hit) # its in the floor
        
        #before_all_hits = True # if all the hits happen after the vertex, veto the event
        
        all_hits_in_floor_wall_before_all_vertices = False
        
        if len(veto_hits) != 0: # there are hits in the floor or wall, consider for veto
            #print('I have veto hits ',len(veto_hits))
            all_hits_in_floor_wall_before_all_vertices = True
        
        if max_vertex_t != -1e6: # a vertex was found
            for hit in veto_hits: # loop over floor and wall hits
                '''
                if self.event_number == 4001:
                    print('vertex time ',max_vertex_t) 
                    print(self.tree.Digi_time[hit])
                '''
                if self.tree.Digi_time[hit] < max_vertex_t:
                    if self.event_number == 4001:
                        print('Saved with ',self.tree.Digi_time[hit])
                    all_hits_in_floor_wall_before_all_vertices = False
                    break
                
                    
        #self.current_vector[self.n] = before_all_hits
        self.current_vector[self.n] = all_hits_in_floor_wall_before_all_vertices
        #self.current_vector[self.n] = 1
    
                
        
        
class pair_track:
    
    def __init__(self, v_hit, t_hit):
        self.hits = [v_hit, t_hit] # v_hit happens before t_hit
        self.xs = [np.zeros(4), np.zeros(4)] # 4 vectors of hits
    
        
    def get_st_points(self, tree):
        # get points in spacetime for each hit
    
        for i, hit in enumerate(self.hits):
            self.xs[i][0] = tree.Digi_time[hit] # i = 0 < i = 1 in time
            self.xs[i][1] = tree.Digi_x[hit]
            self.xs[i][2] = tree.Digi_y[hit]
            self.xs[i][3] = tree.Digi_z[hit]
            
        
    def find_beta(self):
        
        '''
        # we aren't looking for the interval
        ds = (self.xs[0] - self.xs[1])**2 # interval in mostly plus and euclidean time
        
        self.beta = ds[0] / np.sum(ds[1:])
        '''
        self.find_velocity()
        
        self.beta = np.sqrt(np.sum(self.v**2)) / physics.c
        
        return self.beta
        
        
    def find_velocity(self):
        self.x0 = self.xs[0] # initial hit (in floor or wall)
        
        self.v = self.xs[1][1:] - self.x0[1:] # = \Delta \vec{x}
        self.v /= self.xs[1][0] - self.x0[0] # / \Delta t 
        
        
    def find_closest_approach(self, tree, track):
        # see /MATHUSLA-Kalman-Algorithm/docs/Closest_Approach_of_Parametric_Vectors.pdf
        
        v_tr = np.zeros(3)
        x_tr = np.zeros(4)
        
        #x_tr[0] = tree.Track_k_m_t0[track]
        x_tr[1] = tree.Track_k_m_x0[track]
        x_tr[2] = tree.Track_k_m_y0[track]
        x_tr[3] = tree.Track_k_m_z0[track]
        
        v_tr[0] = tree.Track_k_m_velX[track]
        v_tr[1] = tree.Track_k_m_velY[track]
        v_tr[2] = tree.Track_k_m_velZ[track]
        
        del_v = v_tr - self.v
        del_x = x_tr[1:] - self.x0[1:]
        
        t_CA = np.sum(del_x * del_v)
        t_CA /= np.sum(del_v * del_v)
    
        self.closest_approach = np.sqrt(abs(np.sum(del_x + del_v * t_CA))) # square root of the squared distance of closest approach
        
        

class scissors():
    '''since this is what does the cutting!'''

    def __init__(self,file):

        self.file = file
        self.plotter = Plotter() # to be overwritten by passed plotter


    def store_space(self):

        space = sample_space(self.file)
        #space.plotter = self.plotter
        
        self.cut_vectors = space.get_states_m()
        self.func_dicts = space.cuts
        
        #self.plotter = space.plotter # update the plotter to collect information from the parameter space
        #self.plotter.merge_a_plotter(space.plotter) # merge the plotter in space with current plotter to gather data
        self.plotter.event_infos = space.plotter.event_infos
        
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
        
       
        
class event_info:
        ''' An object whose sole purpose is store store information 
            from an event in the analysis. This can be used to access
            information internal to the analysis for plotting. '''
            
        def __init__(self, event):
        
            self.event = event
            
            self.data_dict = {'Exp Pos x wall':[], 'Exp Pos y wall':[],
                              'Exp Pos x floor':[], 'Exp Pos z floor':[],
                              'dt wall':[], 'dt floor':[],
                              'dist wall':[], 'dist floor':[],
                              'vert pos':[], 'slope dist':[],
                              'pair closest approach':[], 'pair reco beta':[],
                              'min hits per track':[], 'max hits per track':[]} # add whatever you want to look at here (and in plotter)!



class Plotter:
    ''' This class organizes the event_info objects into a format usable for plots.
        It's utility is in organizing and accessing event_infos across files. '''
    

    def __init__(self):
        ''' Assign to attributes to collect data for plots '''
        
        self.data_dict = {'Exp Pos x wall':[], 'Exp Pos y wall':[],
                          'Exp Pos x floor':[], 'Exp Pos z floor':[],
                          'dt wall':[], 'dt floor':[],
                          'dist wall':[], 'dist floor':[],
                          'vert pos':[], 'slope dist':[],
                          'pair closest approach':[], 'pair reco beta':[]}
        
        
        self.event_infos = [] # collected info for plotting organized by event
        self.lists_of_event_infos = [] # list of event_infos lists above (one for each file)
        self.survivors = [] # list of survivor inds (one for each file)
        
        self.sample = 'sim data' # type of sample being used to make plots
        
    
    def merge_a_plotter(self, pltr):
       
       #try:
       #for key in pltr.data_dict.keys():
       #    self.data_dict[key].extend(pltr.data_dict[key])
           
       self.lists_of_event_infos.append(pltr.event_infos) # add the event infos as lists with survivors so they can be referenced easily
       self.survivors.append(pltr.survivors)
           
    
    def gather_info(self, key, cut='all'):
        ''' specify which data in event_info.data_dict you want to view
            specifying cut will only look at events at a particular cut,
            otherwise it will gather all entries '''
            
        gathered_info = []
        gathered_event_infos = []
        
        # **** might not need to distinguish because of how we merge now?
        
        #if True:
        try:
            if len(self.lists_of_event_infos) != 0: # passed a whole object (ie. this is a global plotter)
                if type(cut) == int: 
                    for i, event_infos in enumerate(self.lists_of_event_infos):
                        event_infos = np.array(event_infos, dtype=object)
                        survivors = np.array(self.survivors[i][cut], dtype=int)
                          
                        gathered_event_infos.extend(list(event_infos[survivors]))
                        
                elif cut == 'all':
                    gathered_event_infos =  [item for sublist in self.lists_of_event_infos for item in sublist] 
            
            else: # this is for one file (ie. local plotter for one file)      
                if type(cut) == int:
                    event_infos = np.array(self.event_infos, dtype=object)
                    
                    print(self.survivors)
                    
                    survivors = np.array(self.survivors[cut], dtype=int)
                     
                    gathered_event_infos.extend(list(event_infos[survivors]))
                    
                elif cut == 'all':
                    gathered_event_infos = [item for sublist in self.event_infos for item in sublist] 
            
        except IndexError:
            print(event_infos)
            print(survivors)
            print('sorry, I dont think those events have info')
          
        for event_info in gathered_event_infos:
            gathered_info.extend(event_info.data_dict[key])
            
        return gathered_info
    
        
    def plot(self):
        ''' Generate the plots you want with any of the methods! Comment out method calls to ignore the plot '''
    
        #self.plot_min_exp_dist()
        
        #self.expected_hit_positions('wall')
        #self.expected_hit_positions('floor')
        
        #self.expected_time_space('wall')
        #self.expected_time_space('floor')
        
        #self.vertex_pos()
        #self.slope_dist()
    
        self.beta_vs_approach()
        
    
    def beta_vs_approach(self):
        
        _Title = 'Pair Reconstructed Beta vs Closest Approach with Vertex Track {}'.format(self.sample)
        _xlabel = 'closest approach to vertex track [cm]'
        _zlabel = 'pair reconstructed beta []'
        _fname = 'beta_approach_{}.png'.format(self.sample)

        _data_x = np.array(self.gather_info('pair closest approach',2)) #self.data_dict['pair closest approach']
        _data_z = np.array(self.gather_info('pair reco beta',2)) #self.data_dict['pair reco beta']

        #print(_data_x)
        #print(_data_z)
        
        _xlims = [0,np.max(_data_x)*1.1]
        _zlims = [0,1.0] # [np.amin(_data_z),np.max(_data_z)*1.1]
        
        _xbins=100
        _zbins=100
        
        visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, zbins=_zbins,
	      	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
                              
    
    def slope_dist(self):    
        
        _Title = 'Slope Dist ({})'.format(self.sample)
        _xlabel = 'dist [cm]'
        _fname = 'slope_dist_{}.png'.format(self.sample)

        _data_x = self.gather_info('slope dist',-1)
        
        _xlims = [-5000,5000]
        
        _xbins=100
        
        visualization.root_Histogram(_data_x, rng=_xlims, ft_rng=None, bins=0, Title=_Title, xaxis=_xlabel, logx=False, logy=False, fname=_fname)
                   
                
    def plot_min_exp_dist(self):
     
        _Title = 'Minimum Distances to Expected Hit of Digi Hit Before Vertex {}'.format(self.sample)
        _xlabel = 'Minimum Distance to a floor hit [cm]'
        _zlabel = 'minimum distance to a wall hit [cm]'
        _fname = 'min_dist_{}.png'.format(self.sample)

        _data_x = self.data_dict['dist floor']
        _data_z = self.data_dict['dist wall']

        det = detector.Detector()
        _xlims = [np.amin(_data_x),np.max(_data_x)*1.1]
        _zlims = [np.amin(_data_z),np.max(_data_z)*1.1]
        
        _xbins=100
        _zbins=100
        
        visualization.root_Histogram(_data_z, rng=None, ft_rng=None, bins=0, Title=_Title+" (wall)", xaxis=_zlabel, logx=False, logy=False, fname='min_dist_wall_{}.png'.format(self.sample))
        visualization.root_Histogram(_data_x, rng=None, ft_rng=None, bins=0, Title=_Title+" (floor)", xaxis=_xlabel, logx=False, logy=False, fname='min_dist_floor_{}.png'.format(self.sample))
        
        visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, zbins=_zbins,
	      	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
                                                    
        
    def expected_time_space(self, plane):    
        
        _Title = 'Distance to hit vs time to hit in {} ({})'.format(plane, self.sample)
        _xlabel = 'dist [cm]'
        _zlabel = 'time [ns]'
        _fname = 'dist_to_hit_{}_{}.png'.format(self.sample, plane)

        if plane == 'wall':
            _data_x = self.data_dict['dist wall']
            _data_z = self.data_dict['dt wall']
            
        if plane == 'floor':
            _data_x = self.data_dict['dist floor']
            _data_z = self.data_dict['dt floor']
            
        _xlims = [0,2000]
        _zlims = [0,300]
        
        _xbins=100
        _zbins=100
        
        visualization.root_Histogram(_data_z, rng=_xlims, ft_rng=None, bins=0, Title='Time to hit in {} ({})'.format(plane, self.sample), 
                        xaxis=_zlabel, logx=False, logy=True, fname='time_to_hit_{}_{}.png'.format(self.sample,plane))
        
        visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, zbins=_zbins,
	      	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
                   
                              
    def expected_hit_positions(self, plane):
        
        _Title = 'Expected Hit Positions in the {} ({})'.format(plane, self.sample)
        _xlabel = 'x expected positions [cm]'
        _zlabel = '{} expected positions [cm]'.format(['y','z'][plane=='floor'])
        _fname = 'before_dist_{}_{}.png'.format(self.sample, plane)
        
        if plane == 'wall':
            _data_x = self.data_dict['Exp Pos x wall']
            _data_z = self.data_dict['Exp Pos y wall']
            
        if plane == 'floor':
            _data_x = self.data_dict['Exp Pos x floor']
            _data_z = self.data_dict['Exp Pos z floor']
            
        _xlims = [-3000, 5000]
        _zlims = [-3000, 5000]
        _xbins=100
        _zbins=100
            
        visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, zbins=_zbins,
    	        	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
                                                             
            
    def vertex_pos(self):
     
        _Title = 'Fiducial Vertex Positions {}'.format(self.sample)
        _xlabel = 'z position [cm]'
        _zlabel = 'y position [cm]'
        _fname = 'vert_pos_{}.png'.format(self.sample)

        if len(self.data_dict['vert pos']) != 0:
            vertex_pos = np.array(self.data_dict['vert pos'])
        
        else:
            vertex_pos = np.array(self.gather_info('vert pos',1))
        
        joblib.dump(vertex_pos,'vert_pos_array.joblib')

        if vertex_pos.size == 0:
            print('no surviving vertices found to plot!')
            return

        _data_x = vertex_pos[:,2]
        _data_z = vertex_pos[:,1]

        det = detector.Detector()
        _xlims = det.BoxLimits[2]
        _zlims = det.BoxLimits[1]
        
        _xbins=100
        _zbins=100
        
        visualization.root_2D_Histogram(_data_x, _data_z, Title=_Title, xbins=_xbins, zbins=_zbins,
	      	              xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname=_fname)
            
            

def main():
  
    if len(sys.argv) > 1:
        # default set for job submission with files and write locations passed as arguments
        #
        # DO NOT TOUCH!
        
        if len(sys.argv) != 4:
            print('Usage: python cutflow.py <tracker tree> <save directory> <integer identifier>')
            return
    
        plot_cut = False
        plot_obj = False
        
        sum_flows = True # True <=> Background
        load = False
        save = True
        start_from_cut = False
        
        files = [sys.argv[1]]

        #-------------------------

    else:
        # plotting booleans (do you want to make plots?)
        plot_cut = True
        plot_obj = False
        
        sum_flows = True # True <=> Background / sum over data in files for flows ***** need to adress sum_flows or load booleans in below code
        
        load = False # set at most one of these to True
        save = False
        
        start_from_cut = False
        
        start_cut = -1 # work only with the files with survivors at cut start_cut (indexed as in flows)
        # it would be good to have it run without cutting, just to gather passed_event joblibs into a single file
        # add boolean for that so we can start after the relevant cut. 
        
        if start_from_cut:
            #passed_events_file = 'passed_events.joblib'
            passed_events_file = 'passed_events_1_left_28_2_22.joblib'
            #passed_events_file = 'passed_events_run6_4hits_23_2_22.joblib'
            #passed_events_file = 'passed_events_0_left_27_2_22.joblib'
            #passed_events_file = 'passed_events_1_left.joblib'

            passed_events_prev = joblib.load(passed_events_file)
            
            
            # **** could save the scissors after cutting and just update when loaded and cut again
            # **** could also point to the joblib with passed events
            # **** 
            
            # *** use file converter to convert between root and joblib file to start from cut on a load run

            files = []
            
            if load:
                file_converter_file = 'file_converter_W_27_2_22.joblib'
                
                file_converter = joblib.load(file_converter_file) # converts root file to joblib file
            
            for file in passed_events_prev.keys():
                if len(passed_events_prev[file][start_cut]) != 0: # check if it has any surviving events at start_cut
                    if load:
                        files.append(file_converter[file])
                        
                    else:
                        files.append(file)
                        
        elif load:
            if sum_flows:
                #load_files_dir = "/home/keeganh/scratch/job_test/W_sample_dir/run6/analysis_data/21_01_22/"
                #load_files_dir = "/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run3/analysis_data/30_01_22/"
                #load_files_dir = "/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run6/analysis_data/06_02_22/"
                #load_files_dir = "/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run6/analysis_data/12_02_22/"
                #load_files_dir = "/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run6/analysis_data/22_02_22/"
                load_files_dir = "/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run6/analysis_data/27_02_22/"

                files = [filename for filename in glob.iglob(load_files_dir+'/**/scissor_W_*.joblib', recursive=True)]
            
            else:
                load_files_dir = "/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/analysis/save_files/"
                #load_files_dir = "/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/analysis/save_files_run6/"
            
                files = [load_files_dir+'scissor_{}.joblib'.format(sample) for sample in ['h10_0','h2_1','qq_2']]
         
        else:
            if sum_flows: # Background
                directory_3 = '/home/keeganh/scratch/job_test/W_sample_dir/run3/18_12_21/11_24_04/trees/'
                #directory_4 = '/home/keeganh/scratch/job_test/W_sample_dir/run4/09_01_22/09_11_57/trees/'
                #directory_6 = '/home/keeganh/scratch/job_test/W_sample_dir/run6/tracker_data/'
                
                #files = [filename for filename in glob.iglob(directory_6+'/**/stat_*.root', recursive=True)]
                files = [filename for filename in glob.iglob(directory_3+'stat_*.root', recursive=True)][:50]
                #files.extend([filename for filename in glob.iglob(directory_4+'stat_*.root', recursive=True)])
        
            else: # signal
                files_4 = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/21_02_22/13_52_17/trees/stat_0_0.root',
                    '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/21_02_22/13_52_17/trees/stat_1_0.root',
                    '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/21_02_22/13_52_17/trees/stat_2_0.root'] # 4 hits per track
                    
                files_3 = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_0_0.root',
                    '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_2_0.root',
                    '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm//08_01_22/17_50_20/trees/stat_0_0.root'] # 3 hits per track
        #            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_4_0.root']
        
                files = files_4
           
    if len(files) == 0:
        print("I need at least 1 file to run!")
        return

    cuts_to_plot = [9,14] # [int] which cuts to plot (index in cut_options) # **** need to update for multiple in list to work for sum_flows = True as well
    
    sum_values = [[] for i in range(ncuts)] # a list to add values to plot for every cut
    drawers = [] # to hold the scissors

    flows = dict()
    passed_events = dict()
    file_converter = dict()
    
    option = ['cut name', 'cut parameter', 'units', 'on?'] # [name, value of the cut, units of what is being cut on, whether the cut is used (1 => yes)]
    func_name = 'func_name' # func_name is the corresponding data collection function in get_states_m

    # specify cuts and parameters
    cut_options = {'0' :{option[0]:'2 Tracks'                     ,option[1]:2       ,option[2]:'tracks'        ,option[3]:1 , func_name:'Tracks' },
                   '1' :{option[0]:'Vertices'                     ,option[1]:1       ,option[2]:'verts'         ,option[3]:1 , func_name:'Vertices' },
                   '2' :{option[0]:'Fiducial Vertex'              ,option[1]:1       ,option[2]:'bool'          ,option[3]:1 , func_name:'Fiducial_vertex' },
                   '3' :{option[0]:'Floor/Wall Hits Before Vertex',option[1]:2000    ,option[2]:'cm'            ,option[3]:1 , func_name:'Floor_hits_before_vertex' },
                   '4' :{option[0]:'No track hits in the floor'   ,option[1]:True    ,option[2]:'bool'          ,option[3]:1 , func_name:'Track_floor_hits' },
                   '5' :{option[0]:'Vertex opening angle'         ,option[1]:0.02    ,option[2]:'rad'           ,option[3]:0 , func_name:'Opening_angle' },
                   '6' :{option[0]:'Topological Veto'             ,option[1]:-1e2    ,option[2]:'sigma'         ,option[3]:0 , func_name:'Topological' },
                   '7' :{option[0]:'2 Good Betas in a Vertex'     ,option[1]: 1 / 0.2,option[2]:'1 / beta res'  ,option[3]:0 , func_name:'Vertex_track_beta' },
                   '8' :{option[0]:'Hit Differences'              ,option[1]: 1 / -1 ,option[2]:'1 / hits'      ,option[3]:0 , func_name:'Track_hit_diffs' },
                   '9' :{option[0]:'Expected hit edge distance'   ,option[1]:600     ,option[2]:'cm'            ,option[3]:1 , func_name:'Exp_hits' },
                   '10':{option[0]:'No Floor Hits'                ,option[1]:1       ,option[2]:'bool'          ,option[3]:0 , func_name:'No_floor_hits' },
                   '11':{option[0]:'Missing Hit Sum'              ,option[1]:-6      ,option[2]:'-missed hits'  ,option[3]:0 , func_name:'Missing_hit_sum' },
                   '12':{option[0]:'Chi sum cut'                  ,option[1]:-20     ,option[2]:'-chi ndof'     ,option[3]:0 , func_name:'Chi_ndof_cut' },
                   '13':{option[0]:'Hits per track in vertex'     ,option[1]:4       ,option[2]:'hits'          ,option[3]:0 , func_name:'Hits_per_track'},
                   '14':{option[0]:'Fiducial Leniency'            ,option[1]:750     ,option[2]:'cm'            ,option[3]:1 , func_name:'Fiducial_leniency'},
                   '15':{option[0]:'Closest Approach Delta ray'   ,option[1]:750     ,option[2]:'cm'            ,option[3]:1 , func_name:'Delta_ray_cut'},
                   '16':{option[0]:'All vertices before all hits' ,option[1]:1       ,option[2]:'bool'          ,option[3]:0 , func_name:'Fiducial_vertex_time'}} 

    for opt in option:
        flows[opt] = np.array([])

    for cut_ind in cut_options.keys():
        for i in range(len(option)):
            flows[option[i]] = np.append(flows[option[i]], cut_options[cut_ind][option[i]])
            
    flows[option[3]] = flows[option[3]].astype(int)

    permutation = [0,1,2,4,3,15,14,9,13,11,5,6,7,8,10,12,16] # order in which the cuts are performed
                                                    # describes a permutation of 
                                                    # (0, ..., ncuts-1); keys of cut_options
                                                    
    cuts_to_plot = [np.where(np.array(permutation) == cut)[0][0] for cut in cuts_to_plot] # *** find a slicker way to handle cut plotting -> encorporate into plotter?
    
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
            
            print('file: ',i)
        
            scissor = scissors(file)
            
            if plot_obj:
                if not sum_flows:
                    pltr = Plotter()
                
                pltr.sample = sample
                
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
                write_dir = 'save_files'
                #write_dir = 'save_files'
            
                if not os.path.exists(write_dir):
                    os.makedirs(write_dir)
                    
                joblib.dump(scissor, write_dir+'/scissor_{}_{}.joblib'.format(sample,i)) # save the scissor objects
                                                                                                       # => only need to read data once

        if load:
            if i % 100 == 0:
                print(i," files processed")
            
            try:
                scissor = joblib.load(file)
                
            except (FileNotFoundError, EOFError, OSError) as error:
                print("Sorry, I encountered an error with that file")
                i += 1
                continue
                
        scissor.cut_dict(cut_options, permutation)
        
        if not plot_obj: # to save memory (for large runs we run out of memory keeping track of both plot info and survivors)
            passed_events[scissor.file] = scissor.survivor_inds # in order of permutation (indexed as in flows)
            
        if load:
            file_converter[scissor.file] = file # make a file that points between joblib and root file
            file_converter[file] = scissor.file
               
        if load:
            if plot_obj:
                if not sum_flows:
                    pltr = scissor.plotter
            
                pltr.sample = sample
                
                scissor.plotter.survivors = scissor.survivor_inds
                
                pltr.merge_a_plotter(scissor.plotter)
                
                if not sum_flows:
                    pltr.plot()
    
        else:
            if plot_obj:
                scissor.plotter.survivors = scissor.survivor_inds

                # *** right now it's storing survivor inds as a duplicate; we should find a way that we don't have to do this;
                # ie. just store them in the plotter and not in the main scissor object?
                
                pltr.merge_a_plotter(scissor.plotter)
            
            if plot_obj and not sum_flows:
                pltr.plot() # generate plots for this sample with plotter
            
    
        file_flows = scissor.flows.astype(int)
        file_flows = file_flows.take(permutation,0) # permute the flows to the order the cuts were performed
        
        if sum_flows: # and not plot_obj:  *** implement a seperate boolean for deciding whether or not to show flows.  (leave off for job default)
            try:
                flows['Flows'] += file_flows
                             
            except:
                flows['Flows'] = file_flows

       
        else:
            flows[sample] = file_flows
            flows[sample+' fraction'] = np.round(file_flows / scissor.events, 4)
        
        total_events += scissor.events
 

        if plot_cut and len(scissor.survivor_inds) != 0:
            #------- Plot the Chosen cuts
            for cut in cuts_to_plot:
                inds = np.array(scissor.survivor_inds[cut-1],dtype=int) # surviving indices before cut
                
                values = scissor.cut_vectors[inds,scissor.func_dicts[cut_options[str(permutation[cut])]['func_name']]['index']]
    
                if sum_flows:
                    sum_values[cut].extend(values) # ****** need to change this when we are plotting more than one thing!!! 
    
                if not sum_flows and len(values) != 0:
                    _bins = 100
                    #_rng = (np.amin(values),np.max(values)*1.1)
                    _rng = (0,5000)
                    visualization.root_Histogram(values,
                  					rng=_rng,
                  					bins=_bins-int(np.sqrt(len(values))),
                  					Title='{} {}'.format(flows['cut name'][cut],sample),
                  					xaxis=flows['units'][cut],
                  					fname='distribution_{}_{}.png'.format(cut,sample))
                    
            #-----------
                
        i += 1

    if sum_flows:
        
        if not plot_obj: # to save memory
            flows['Flows fraction'] = np.round(flows['Flows'] / total_events, 8)
        
        print("Total Events: ",total_events)

        if plot_obj:
            pltr.plot() # generate plots for the run with plotter
        
        if plot_cut and len(sum_values) != 0:
            for cut in cuts_to_plot:
                _bins = 100
                _rng = (np.amin(sum_values[cut]),np.max(sum_values[cut])*1.1)
                if cut == 4:
                    _rng = (-3000,5000)
                elif cut == 7:
                    _rng = (0,3000)
                visualization.root_Histogram(sum_values[cut],
              					rng=_rng,
              					bins=_bins-int(np.sqrt(len(values))),
              					Title='{} {}'.format(flows['cut name'][cut],sample),
              					xaxis=flows['units'][cut],
              					fname='distribution_{}_{}.png'.format(cut,sample),
                        logy=False,
                        logx=False)
                    
    cutflow = pd.DataFrame(flows)

    print(cutflow)
    
    if len(sys.argv) == 1: # and not start_from_cut:
        joblib.dump(passed_events,'passed_events.joblib')
    
    if load:
        joblib.dump(file_converter,'file_converter.joblib')

    if plot_obj:
        joblib.dump(pltr,'global_plotter.joblib')

    
    
if __name__ == '__main__':

    main()



