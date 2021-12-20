#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import ROOT as root
import numpy as np


for n in range(3,11):
	found = False

	for chi in np.linspace(70,200,50):
		ndof = 4 * n - 6

		if root.Math.chisquared_cdf(chi,ndof) == 1 and not found:
			print("{} ndof has chi threshold {} and chi ndof threshold {}".format(ndof,chi,chi/ndof))

			found = True
