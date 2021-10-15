#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 19:11:39 2021

@author: keegan
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np


L = 3 # [cm]
L_Al = 1 # [cm]
L_Sc = 2 # [cm]
L_r_Al = 24.0111 # [g cm^(-2)] Aluminum
L_r_Sc = 43 # [cm] Scintillator (Saint-Gobain paper)

#p = 4.85e4 # [MeV] from muon momentum distribution
p = 3e3


L_r_lyn = L_Al / L_r_Al + L_Sc / L_r_Sc

sigma_lyn = 13.5 * np.sqrt(L_r_lyn) * (1 + 0.038 * np.log(L_r_lyn))

print("M.S. STD [rad MeV] Lynch Prescription is: {:.03g} rad".format(sigma_lyn)) 

print("M.S. STD [rad] is: {:.03g} rad".format(sigma / p)) 

