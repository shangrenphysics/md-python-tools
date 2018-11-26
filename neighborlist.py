#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 13:38:17 2018

@author: shangren

This script is designed to establish the neighbor list within the cutoff
Input parameter sigma is the cutoff of neighbor list
Input file is any dump file of lammps. Default form of coordinates is scaled, i.e. xs, ys, etc

"""

import numpy as np
import sys


def neighborlist(ai,sigma,fin):
   data = np.loadtxt(fin,skiprows=9)
   a = data[ai-1][0:5]
   shift = np.subtract(a[2:5],[0,0,0])
   cor_list = np.subtract(data[:,2:5],shift)
   cor_list[cor_list<-0.5] = cor_list[cor_list<-0.5]+1
   cor_list[cor_list>0.5] = cor_list[cor_list>0.5]-1
   distance=np.sqrt(np.diagonal(np.dot(cor_list,np.transpose(cor_list))))
   indices_sorted=np.argsort(distance)
   output=np.zeros((np.shape(data)[0],6))
   output[:,0:5] = data[indices_sorted,0:5]
   output[:,5] = distance[indices_sorted]
   return output

fin = sys.argv[1]
print(neighborlist(2,1,fin))
