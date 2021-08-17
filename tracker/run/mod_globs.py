#!/usr/bin/env python

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import numpy as np
from joblib import dump, load



class Editor:

        def __init__(self, filename, i, j):

                self.pars = np.genfromtxt('Parameters.txt',delimiter=',')

                if len(np.shape(self.pars)) == 1: #1 dimensional array not two
                        self.pars = [self.pars]

                self.i = i
                self.j = j

                self.fname = filename
                self.f = open(self.fname,"r+")

                self.lines = self.f.readlines()

        def alter(self):

                j = 0
                for line in self.lines:
                        if "kalman_chi_s" in line:
                                self.lines[j] = "        const double kalman_chi_s = {};\n".format(self.pars[self.i][0])

                        elif "kalman_chi_add" in line:
                                self.lines[j] = "        const double kalman_chi_add = {};\n".format(self.pars[self.i][1])

                        elif "kalman_track_chi" in line:
                                self.lines[j] = "        const double kalman_track_chi = {};\n".format(self.pars[self.i][2])

                        elif "double p " in line:
                                self.lines[j] = "        const double p = {}; // [MeV] representative momentum\n".format(self.pars[self.i][3])

                        elif "merge_cos_theta" in line:
                                self.lines[j] = "        const double merge_cos_theta = {};\n".format(self.pars[self.i][4])

                        elif "merge_distance" in line:
                                self.lines[j] = "        const double merge_distance = {}*units::cm;\n".format(self.pars[self.i][5])

                        elif "seed_closest_approach" in line:
                                self.lines[j] = "        const double seed_closest_approach = {}*units::cm;\n".format(self.pars[self.i][6])

                        elif "vertex_chi2" in line:
                                self.lines[j] = "        const double vertex_chi2 = {};\n".format(self.pars[self.i][7])

                        elif "closest_approach_add" in line:
                                self.lines[j] = "        const double closest_approach_add = {}*units::cm;\n".format(self.pars[self.i][8])

                        elif "kalman_vertex_chi_add" in line:
                                self.lines[j] = "        const double kalman_vertex_chi_add = {};\n".format(self.pars[self.i][9])

                        elif "kalman_vertex_chi" in line:
                                self.lines[j] = "        const double kalman_vertex_chi = {};\n".format(self.pars[self.i][10])

                        elif "kalman_v_add" in line:
                                self.lines[j] = "        const std::vector<double> kalman_v_add = {%f,%f};\n"%(self.pars[self.i][11],self.pars[self.i][12])

                        elif "kalman_v_drop" in line:
                                self.lines[j] = "        const std::vector<double> kalman_v_drop = {%f,%f};\n"%(self.pars[self.i][13],self.pars[self.i][14])

                        elif "start_ev" in line:
                                self.lines[j] = "        const int start_ev = {:.00f};\n".format(self.pars[self.i][15 + 2 * self.j])

                        elif "end_ev" in line:
                                self.lines[j] = "        const int end_ev = {:.00f};\n".format(self.pars[self.i][16 + 2 * self.j])


                        j += 1

        def write(self):

                self.f.close()
                self.f = open(self.fname,"w")

                self.f.writelines(self.lines)



def main():

	i = int(sys.argv[1])
	j = int(sys.argv[2])

	print("current parameter counter: ",i)
	print("current dataset counter: ",j)

	edt = Editor("../include/globals.hh",i,j)

	print(edt.pars[i])

	edt.alter()
	edt.write()

	if i == np.shape(edt.pars)[0]-1: #last parameter has been used, increment j

		g = open("done.txt","w+")
		g.write("All done! Please delete me.")
		g.close()



if __name__ == "__main__":

	main()

