#!/usr/bin/env python

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import numpy as np
from joblib import dump, load



class Editor:

        def __init__(self, filename, i, j):

                self.pars = np.genfromtxt('Parameters.txt',delimiter=',')

                if len(np.shape(self.pars)) == 1: #1 dimensional array, need 2
                        self.pars = [self.pars]

                self.i = i
                self.j = j

                self.fname = filename
                self.f = open(self.fname,"r+")

                self.lines = self.f.readlines()

        def alter(self):

                pars_to_handle = ['kalman_chi_s',
				'kalman_chi_add',
				'kalman_track_chi',
				'kalman_pval_drop',
				'kalman_pval_add',
				'kalman_pval_track',
				'p',
				'merge_cos_theta',
				'merge_distance',
				'seed_closest_approach',
				'vertex_chi2',
				'closest_approach_add',
				'kalman_vertex_chi_add',
				'kalman_vertex_chi',
				'kalman_v_add[0]',
				'kalman_v_add[1]',
				'kalman_v_drop[0]',
				'kalman_v_drop[1]',
				'start_ev',
				'end_ev']

                pars_dict = dict()

                k = 0
                for par in pars_to_handle:
                        pars_dict[par] = self.pars[self.i][k]
                        k += 1

		# start and end event for each dataset
                #pars_dict[pars_to_handle[-2]] = self.pars[self.i][15 + 2 * self.j]
                #pars_dict[pars_to_handle[-1]] = self.pars[self.i][16 + 2 * self.j]
                pars_dict[pars_to_handle[-2]] = self.pars[self.i][-2]
                pars_dict[pars_to_handle[-1]] = self.pars[self.i][-1]

                for k in range(len(self.lines)):
                    par = self.lines[k].split(' ')[0]

                    if par in pars_to_handle:
                        self.lines[k] = par + ' {}\n'.format(pars_dict[par])

        def write(self):

                self.f.close()
                self.f = open(self.fname,"w")

                self.f.writelines(self.lines)



def main():

	i = int(sys.argv[1])
	j = int(sys.argv[2])

	print("current parameter counter: ",i)
	print("current dataset counter: ",j)

#	edt = Editor("../include/globals.hh",i,j)
	edt = Editor("par_card.txt",i,j)

#	print(edt.pars[i])

	edt.alter()
	edt.write()

	if i == np.shape(edt.pars)[0]-1: #last parameter has been used, increment j
		g = open("done.txt","w+")
		g.write("All done! Please delete me.")
		g.close()



if __name__ == "__main__":

	main()

