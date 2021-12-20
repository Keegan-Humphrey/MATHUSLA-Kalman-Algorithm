
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

ncuts = 11



class sample_space():


    def __init__(self,file):

        self.file = file
        self.det = detector.Detector()

    def inside_box(self, x, y, z):

        # box_lims = [[-5000., 5000.], [5900.0, 9000.0], [6900., 17100.]]
        #box_lims = [[-5000., 5000.], [6000.0, 8917.0], [7000., 17000.]]
#        box_lims = detector.BoxLimits
        box_lims = self.det.BoxLimits
        y_floor = self.det.LayerYLims[2][1]

        if box_lims[0][0] < x and box_lims[0][1] > x:
           if y_floor < y and box_lims[1][1] > y: # only accept it if above the floor layers
#            if box_lims[1][0] < y and box_lims[1][1] > y:
                if box_lims[2][0] < z and box_lims[2][1] > z:
                    return True

        return False


    def in_layer(self,y_val):

        #LAYERS_Y=[[6001.0, 6004.0],  [6104.0, 6107.0], [8001.0, 8004.0], [8104.0, 8107.0], \
        #[8501.0, 8504.0], [8604.0, 8607.0], [8707.0, 8710.0], [8810.0, 8813.0], [8913.0, 8916.0]]
        '''
        LAYERS_Y=[[6003.0 + 547, 6006.0 + 547],  #layer 0 (floor)
		[6106.0 + 547, 6109.0 + 547], #layer 1 (floor)
		[6209.0 + 547, 6212.0 + 547], #layer 2 (floor)
 		[8003.0 + 547, 8006.0 + 547], #layer 3
		[8106.0 + 547, 8109.0 + 547], #layer 4
		[8503.0 + 547, 8506.0 + 547], #layer 5
		[8606.0 + 547, 8609.0 + 547], #layer 6
		[8709.0 + 547, 8712.0 + 547], #layer 7
		[8812.0 + 547, 8815.0 + 547], #layer 8
		[8915.0 + 547, 8918.0 + 547]] #layer 9
        '''
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

        self.tree.SetBranchStatus("Digi_y", 1)
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


        cut_vectors = []
        
        #print(self.tree.Vertex_k_m_trackIndices)
#        self.fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)

        self.expected_x_total = np.zeros((int(self.tree.GetEntries()),2))
        self.expected_z_total = np.zeros((int(self.tree.GetEntries()),2))

        for event_number in range(int(self.tree.GetEntries())):

            self.event_number = event_number

            self.current_vector = np.zeros(1+ncuts) # first will be the event number
            self.n = 0

            self.tree.GetEntry(event_number)
            if event_number % 1000 == 0:
                print("event:", event_number)

            # record the event number for reference
            self.current_vector[0] = event_number
            self.n += 1

            #------------- cut on the number of made tracks in the event
            self.current_vector[1] = self.tree.NumTracks_k_m
            self.n += 1

            #------------- cut on the number of vertices made in the event
            self.current_vector[2] = len(self.tree.Vertex_k_m_x)
            self.n += 1
            
            #------------- Optional Cuts
            
            self.Fiducial_vertex()
            self.n += 1
            
            self.Floor_hits_before_vertex()
            self.n += 1

#            self.Expected_hits_old()
            self.n += 1 # add this for every commented out line
            
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

            #------------ All done
            cut_vectors.append(self.current_vector)

        self.cut_vectors = np.array(cut_vectors)

        return self.cut_vectors


    def Fiducial_vertex(self):
            #------------- Fiducial vertex cut
            inside = False

            for num in range(len(self.tree.Vertex_k_m_x)):
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]
                if self.inside_box(vtxx, vtxy, vtxz):
                    inside = True

            # check if at least one vertex is in the detector
            self.current_vector[3] = int(inside) # fiducial
            

    def Floor_hits_before_vertex(self):
            #------------- Floor hits before vertex
            floorveto = False

            '''
            min_digi_t = 1e6
            for hit in range(len(self.tree.Digi_y)): # loop over the hits
                if self.in_layer(self.tree.Digi_y[hit]) <= 2: # check if layer index is a floor layer index
                    if self.tree.Digi_time[hit] < min_digi_t: # find minimum t in the floor
                        min_digi_t = self.tree.Digi_time[hit]

            vetoed_vertices = 0

            for num in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]

                if self.tree.Vertex_k_m_t[num] > self.tree.Digi_time[hit]: # see if any floor hits happened before the vertex
                    vetoed_vertices += 1

                elif not self.inside_box(vtxx, vtxy, vtxz): # if they did, make sure the vertex is fiducial (above the floor)
                    vetoed_vertices += 1

            if len(self.tree.Vertex_k_m_x) == vetoed_vertices:
                floorveto = True

            '''
            for num in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vertexveto = False
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]

                for hit in range(len(self.tree.Digi_y)): # loop over the hits
                    if self.in_layer(self.tree.Digi_y[hit]) <= 2: # check if layer index is a floor layer index
                        if self.tree.Vertex_k_m_t[num] > self.tree.Digi_time[hit]: # check if the floor hit happened before the vertex
                            vertexveto = True

                if not self.inside_box(vtxx, vtxy, vtxz): # make sure the vertex is fiducial (above the floor)
                    vertexveto = True

                if not vertexveto: # There's a vertex in the event that isn't vetoed
                    break

            if len(self.tree.Vertex_k_m_x) != 0:
                if vertexveto:
                    floorveto = True

            # track floor veto (this is really a timing veto)
#            self.current_vector[4] = int(not floorveto) # cut if floorveto == True => 0 == current_vector[4] (< 1 is true)
            self.current_vector[self.n] = int(not floorveto) # cut if floorveto == True => 0 == current_vector[4] (< 1 is true)


    def Expected_hits_old(self):
            #---------- Expected Hits (old)
            expectedveto = False
#            self.fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
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

#            self.current_vector[5] = int(expectedveto) # bottom layer hits
            self.current_vector[self.n] = int(expectedveto) # bottom layer hits


    def Opening_angle(self):
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
			# X cut on smallest opening angle
    #            	self.current_vector[6] = open_angles[0]
                	self.current_vector[self.n] = open_angles[0]
                                                   

    def Topological(self):
            #---------------- check if any vertices are below any track hits within some sigma (Topological Cut)
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

#                pulls.append((self.tree.Vertex_k_m_y[k1] - min_y) / self.tree.Vertex_k_m_ErrorY[k1]) # pull of lowest hit from vertex

                if self.tree.Vertex_k_m_ErrorY[k1] != 0:
                    pulls.append((min_y - self.tree.Vertex_k_m_y[k1]) / self.tree.Vertex_k_m_ErrorY[k1]) # pull of lowest hit from vertex
                else:
                    pulls.append(1e6)

            if len(pulls) != 0:
#                current_vector[7] = np.max(pulls)
                self.current_vector[self.n] = np.max(pulls)


            '''
                        if self.tree.Digi_y[hit] < self.tree.Vertex_k_m_y[k1]:
                            vertexveto = True
                if vertexveto:
                    topoveto += 1

            if topoveto < current_vector[2]: # vetoed vertices < total vertices
                # see if there is at least one vertex that isn't vetoed
                current_vector[7] = 1
            '''


    def Vertex_track_beta(self):
            #---------- Vertex Beta Cut
            # want at least one vertex with two betas within some threshold in each event

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
#                current_vector[9] = 1/np.amin(min_diffs) if np.amin(min_diffs) != 0 else 1e6
                self.current_vector[self.n] = 1/np.amin(min_diffs) if np.amin(min_diffs) != 0 else 1e6


    def Exp_hits(self):
            #------------- Expected Hits Cut

            max_centre_dist = 0

            centre = [np.mean(self.det.BoxLimits[0]), np.mean(self.det.BoxLimits[2])]

            floor_ys = [np.sum(self.det.LayerYLims[i]) / 2 for i in range(3)]
            floor_y = np.sum(self.det.LayerYLims[2]) / 2

            expected_x = []
            expected_z = []

#        floor_y = np.sum(Detector.LayerYLims[2]) / 2

            vertex_trackIndices = self.tree.Vertex_k_m_trackIndices
            vertex_trackIndices = util.unzip(vertex_trackIndices)

            hits_in_floor = False

            for i in range(3):
                if np.any(self.tree.Digi_y == np.sum(self.det.LayerYLims[i]) / 2):
                    hits_in_floor = True

            if len(vertex_trackIndices) > 0 and not hits_in_floor:
                for vert in range(len(vertex_trackIndices)): # vertex
                    vert = int(vert)

                    for ind in vertex_trackIndices[vert]: # track
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


#            self.expected_x_total.extend(expected_x)
#            self.expected_z_total.extend(expected_z)

#        expected_x = np.abs(expected_x)
#        expected_x = np.abs(expected_z)

            edge_dist_x = (self.det.BoxLimits[0][1] - self.det.BoxLimits[0][0]) / 2 - np.abs(np.array(expected_x) - centre[0]) # distance from expected hits
            edge_dist_z = (self.det.BoxLimits[2][1] - self.det.BoxLimits[2][0]) / 2 - np.abs(np.array(expected_z) - centre[1]) # to edge of detector

            if len(edge_dist_x) != 0:
                 min_edge_dist = np.amin(np.amin([edge_dist_x, edge_dist_z], axis=1))

                 self.expected_x_total[self.event_number] = expected_x[:2] #edge_dist_x[:2]
                 self.expected_z_total[self.event_number] = expected_z[:2] #edge_dist_z[:2]

#                 current_vector[10] = min_edge_dist
                 self.current_vector[self.n] = min_edge_dist

            # need a wall veto cut as well


    def No_floor_hits(self):
        ''' Veto the event if it has floor hits '''
    
        floor_veto = False
    
        for hit in range(len(self.tree.Digi_y)): # loop over the hits
            if self.in_layer(self.tree.Digi_y[hit]) <= 2: # check if layer index is a floor layer index
                floor_veto = True
                break

        if floor_veto:
            self.current_vector[self.n] = True  
        

class scissors():
    '''since this is what does the cutting!'''

    def __init__(self,file):

        self.file = file

    def store_space(self):

        self.space = sample_space(self.file)

        # self.cut_vectors = self.space.get_states()
        self.cut_vectors = self.space.get_states_m()

        self.events = len(self.cut_vectors)

        print('number of events is: ',len(self.cut_vectors))


    def cut(self, *p, on_off): #p0, p1, p2, p3, p4, p5, p6, p7, p8):

        global ncuts

        self.cuts = [*p]  #p0, p1, p2, p3, p4, p5, p6, p7, p8]

        self.flows = np.zeros(ncuts)
        local_vectors = np.copy(self.cut_vectors)
        self.survivor_inds = []

        #print(np.shape(local_vectors))

        if local_vectors.size != 0:
            for i in range(ncuts):

                #print('local vectors: ',len(local_vectors))

                if on_off[i] == 1:
                    inds = np.where(local_vectors[:,i+1] < self.cuts[i]) # +1 since first is ev number

                else:
                    inds = np.where(local_vectors[:,i+1] < -1e8) # +1 since first is ev number

                ## just revert to previous by having large negative cut 
                local_vectors = np.delete(local_vectors, inds, 0)

                self.flows[i] += np.shape(local_vectors)[0]
                self.survivor_inds.append(local_vectors[:,0])

        self.survivors = self.flows[-1]
        
        if len(self.survivor_inds) != 0:
            print('Events surviving cuts: ',self.survivor_inds[-1])

#        print("flows are ",self.flows)
#        print("with parameters ",self.cuts)



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



def main():

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




if __name__ == '__main__':

    '''
    main()

    files = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_0_0.root',
            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_2_0.root',
            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/25_11_21/20_54_21/trees/stat_4_0.root']
    '''

#    directory = sys.argv[1]
#    directory = '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/09_12_21/18_42_42/trees/'

#    directory = '/home/keeganh/scratch/job_test/W_sample_dir/run3/18_12_21/11_24_04/trees/'
#    files = [filename for filename in glob.iglob(directory+'stat_*.root', recursive=True)]
    '''
    '''

 
    files = ['/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/13_10_21/08_16_40/trees/stat_0_0.root',
            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/13_10_21/08_16_40/trees/stat_2_0.root',
            '/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/13_10_21/13_20_32/trees/stat_4_0.root']
    '''
	    '/home/keeganh/scratch/presentation_dir/sim_trees/h10_1mill.root',
	    '/home/keeganh/scratch/presentation_dir/sim_trees/h2_1mill.root',
	    '/home/keeganh/scratch/presentation_dir/sim_trees/W_1mill.root']

    '''
#    scissor = scissors(sys.argv[1])

#    files = [sys.argv[1]]

    plot = True
    sum_flows = False
    
    sum_values = []
    sum_data_x = []
    sum_data_z = []

    flows = dict()

    cut_name = ['2 Tracks',
			'Vertices',
			'Fiducial Vertex',
			'Floor Hits Before Vertex',
			'Expected Hits',
			'Vertex opening angle',
			'Topological Veto',
			'2nd Good Betas in a Vertex',
			'Hit Differences',
                        'Expected hit edge distance',
                        'No Floor Hits']
                        
    flows['cut name'] = cut_name

    cut_pars = [2,1,1,1,-1,0.02,-1e6,0.2,-1,100,1] # cut the event if it's cut parameter is less than corresponding cut here
#    cut_units = ['tracks','verts','bool','bool','bool','rad','bool','beta res','hits','cm']
    cut_units = ['tracks','verts','bool','bool','bool','rad','sigma','beta res','hits','cm','bool']

    cut_on_off = [1,1,1,1,0,1,1,1,0,1,0]

    flows['cut parameter'] = np.copy(cut_pars)
    flows['units'] = cut_units
    flows['on?'] = cut_on_off

    cut_pars[7] = 1 / cut_pars[7] # invert if you want greater than instead of less than cut
    cut_pars[8] = 1 / cut_pars[8]

    total_events = 0

    i = 0

    for file in files:


        if sum_flows:
            sample = 'W'

        else:
            sample = ['h10','h2','W'][i]

        print(file)

        scissor = scissors(file)

        scissor.store_space()

#        scissor.cut(2,1,1,1,-1,0.02,-1,1/0.2,1/3)
        scissor.cut(*cut_pars,on_off=cut_on_off)

#        print("flows are ",scissor.flows)
#        print("flow fractions are ",scissor.flows / scissor.events)

        if sum_flows:
            try:
                flows['Flows'] += scissor.flows.astype(int)
                total_events += scissor.events
#                flows['Flows fraction'] += np.round(scissor.flows / scissor.events, 4)
            except:
                flows['Flows'] = scissor.flows.astype(int)
                total_events += scissor.events
#                flows['Flows fraction'] = np.round(scissor.flows / scissor.events, 4)

        else:
            flows[sample] = scissor.flows.astype(int)
            flows[sample+' fraction'] = np.round(scissor.flows / scissor.events, 4)


#        print("events with vertices ",scissor.survivor_inds[1])
#        if i == 2:
#            print("Background survivor inds are ",np.array(scissor.survivor_inds[-1]))

#def plot(scissor, values, cut_name, cut_units, sample='data'):

        if plot and len(scissor.survivor_inds) != 0:

            cut = 6 # which cut to view

            inds = np.array(scissor.survivor_inds[cut-1],dtype=int) # indices before cut

#            values = 1/scissor.cut_vectors[inds,cut+1]

            values = scissor.cut_vectors[inds,cut+1]

            if sum_flows:
                sum_values.append(values)
#            print(values)

            _bins = 100
            _rng = (-1000,1000)

            if not sum_flows:
                visualization.root_Histogram(values,
					rng=_rng,
					#rng=(np.amin(values),np.max(values)*1.1),
					bins=_bins-int(np.sqrt(len(values))),
					Title='{} {}'.format(cut_name[cut],sample),
#					xaxis='minimum edge distance for expected hits (no floor hits) [cm]',
					xaxis=cut_units[cut],
					fname='distribution_{}_{}.png'.format(cut,sample))


            '''
            '''

            _data_x = scissor.space.expected_x_total[inds]
            _data_z = scissor.space.expected_z_total[inds]


#            _data_z = scissor.space.expected_z_total[inds]

            _data_x = _data_x.flatten()
            _data_z = _data_z.flatten()

#            _xlims =  [np.amin(_data_x),np.max(_data_x)]
#            _zlims =  [np.amin(_data_z),np.max(_data_z)]
#            _xlims =  [np.amin(_data_x),np.max(_data_x)]
#            _zlims =  [np.amin(_data_z),np.max(_data_z)]

            _xlabel = 'x displacement from edge [cm]'
            _zlabel = 'z displacement from edge [cm]'
#            _xlabel = 'track 1 z displacement from edge [cm]'
#            _zlabel = 'track 2 z displacement from edge [cm]'


            det = detector.Detector()
            _xlims = [det.BoxLimits[0][0]-3000, det.BoxLimits[0][1]+3000]
            _zlims = [det.BoxLimits[2][0]-3000, det.BoxLimits[2][1]+3000]
#            _xlims = [-1000, (det.BoxLimits[0][1]-det.BoxLimits[0][0])/2]
#            _xlims = [0, (det.BoxLimits[2][1]-det.BoxLimits[2][0])/2]
#            _zlims = [-1000, (det.BoxLimits[2][1]-det.BoxLimits[2][0])/2]

            if not sum_flows:
                visualization.root_2D_Histogram(_data_x, _data_z, Title='Expected Positions (No Hits) {}'.format(sample), xbins=100, zbins=100,
    	        	xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname='expected_2D_{}.png'.format(sample))


#        if plot and not sum_flows:
#            plot(scissor, sample)

        i += 1

    if sum_flows:
        flows['Flows fraction'] = np.round(flows['Flows'] / total_events, 8)

        if plot:
            visualization.root_2D_Histogram(_data_x, _data_z, Title='Expected Positions (No Hits) {}'.format(sample), xbins=100, zbins=100,
    	        	xlims=_xlims, zlims=_zlims, xlabel=_xlabel, zlabel=_zlabel, fname='expected_2D_{}.png'.format(sample))

#        plot(scissor, values, cut_name, cut_units, sample='W')

    import joblib
    
    joblib.dump(flows,'flows.pkl')

    cutflow = pd.DataFrame(flows) #, index=cut_name)

    print(cutflow)

