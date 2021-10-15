\import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import matplotlib.pyplot as plt
import ROOT as root
import sys



def find_missing():

    events = [33]
    layers_to_check = [8004.5, 8107.5]

    files = []

    file_dir = sys.argv[1]

    for file in os.listdir(file_dir):
        if file.endswith(".root"):
            files.append(file_dir + file)

    detector_x_limits = [[-4955, -4045], [-3955, -3045], [-2955, -2045], [-1955, -1045], [-955, -45], [45, 955], [1045, 1955], [2045, 2955], [3045, 3955], [4045, 4955]]
    detector_z_limits = [[6995, 7905], [7995, 8905], [8995, 9905], [9995, 10905], [10995, 11905], [11995, 12905], [12995, 13905], [13995, 14905], [14995, 15905], [15995, 16905]]

    missing_z = []
    missing_y = []
    missing_x = []

    digi_y_list = []

    count = 0.0
    track_count = 0.0

    for file in files:

        tfile = root.TFile.Open(file)
        tree = tfile.Get("integral_tree")

#        for n in range(tree.GetEntries()):
        for n in events:
            tree.GetEntry(n)

            print("Truth location ",tree.GenParticle_y[0]/10,", ",tree.GenParticle_x[0]/10,", ",tree.GenParticle_z[0]/10)

            #if tree.NumTracks == 1:

            digi_y_list.clear()
            for j, digi_y in enumerate(tree.Digi_y):
                digi_y_list.append(digi_y)

            for k, value in enumerate(tree.Track_k_m_velY):
                track_count += 1.0

                    #if value > 0:
                for layer_y in layers_to_check:
                    if layer_y > tree.GenParticle_x[0]/10:
                        #x, y, z, t = tree.Track_x0[k], tree.Track_y0[k], tree.Track_z0[k], tree.Track_t0[k]
                        #px, py, pz = tree.Track_velX[k], tree.Track_velY[k], tree.Track_velZ[k]
                        x, y, z, t = tree.Track_k_m_x0[k], tree.Track_k_m_y0[k], tree.Track_k_m_z0[k], tree.Track_k_m_t0[k]
                        px, py, pz = tree.Track_k_m_velX[k], tree.Track_k_m_velY[k], tree.Track_k_m_velZ[k]

                        first_layer_t = (layer_y - y) / py
                        first_layer_x = x + (first_layer_t*px)
                        first_layer_y = layer_y;
                        first_layer_z = z + (first_layer_t*pz)

                        '''
                        second_layer_t = (6105.5 - y) / py
                        second_layer_x = x + (second_layer_t*px)
                        second_layer_y = 6105.5;
                        second_layer_z = z + (second_layer_t*pz)
                        '''
                        for x_limit in detector_x_limits:
                            for z_limit in detector_z_limits:
                                if first_layer_x < x_limit[1] and first_layer_x > x_limit[0] and first_layer_z < z_limit[1] and first_layer_z > z_limit[0]:
                                    # second_layer_x < x_limit[1] and second_layer_x > x_limit[0] and second_layer_z < z_limit[1] and second_layer_z > z_limit[0]:
                                    if digi_y_list.count(layer_y) == 0: # and digi_y_list.count(6105.5) == 0:
                                        missing_z.append(first_layer_z)
                                        missing_x.append(first_layer_x)
                                        missing_y.append(first_layer_y)
                                        count += 1.0


            print ("Missed count: ", count)
            #print ("Hit count: ", hit)
            print ("Track_Count: ", track_count)
            print ("missing x: ",missing_x)
            print ("missing y: ",missing_y)
            print ("missing z: ",missing_z)
    #plt.scatter(missing_z, missing_x, marker='o', s=0.4, color='black')
    #plt.show()




if __name__ == "__main__":

    find_missing()
