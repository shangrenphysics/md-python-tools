#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:29:17 2018

This script is designed to dig the first nth nearest neigbors by the given neighbor number n
Input parameter n is the # of neighbor atoms
Input file is any dump file of lammps. Default form of coordinates is scaled, i.e. xs, ys, etc

@author: shangren

"""

import numpy as np
import sys


def neighborlist(ai,n,fin): #ai=atom_index sigma=cutoff n=number_of_nearest_neighbor
   data = np.loadtxt(fin,skiprows=9)
   a = data[ai-1][0:5]
   shift = np.subtract(a[2:5],[0,0,0])
   cor_list = np.subtract(data[:,2:5],shift)
   #scale coordinates#
   cor_list[cor_list<-0.5] = cor_list[cor_list<-0.5]+1
   cor_list[cor_list>0.5] = cor_list[cor_list>0.5]-1
   distance=np.sqrt(np.diagonal(np.dot(cor_list,np.transpose(cor_list))))
   indices_sorted=np.argsort(distance)
   output=np.zeros((n,6))
   output[:,0:2] = data[indices_sorted,0:2][0:n]
   output[:,2:5] = cor_list[indices_sorted,0:3][0:n]*Lx
   output[:,5] = distance[indices_sorted][0:n]*Lx
   #print(output)
   return output

def nbl2xyz(nbl,fout1,fout2): # this is a function to convert neighbor list to .xyz file, which can be read by VESTA
    if nbl[0,1] == 1:
        fout = fout1
    if nbl[0,1] == 2:
        fout = fout2
    print(len(nbl),file=fout)
    print(int(nbl[0,0]),file=fout)
    for line in nbl:
        if line[1] == 1:
            aty = 'Al' #the default atom type is 'Al' and 'B'. Change it if needed.
        if line[1] == 2:
            aty = 'B'
        print(aty,*line[2:5],file=fout)
    return 0

def getLx(fin): #Lx is needed to convert the coordinates of cluster from scale format to real format. In our case simulation box is cubic, i.e. Lx=Ly=Lz
    pf = open(fin,'r')
    lines = pf.readlines()
    pf.close()
    global Lx
    Lx = float(lines[5].split()[1])-float(lines[5].split()[0])
    na = int(lines[3].split()[0])
    return Lx, na

#fin = sys.argv[1]
fin = 'all.6996000.atom'
n = 16 # a parameter which means the number of neighbors. Change it according to the coordination number analysis.
Lx, na = getLx(fin)
fout1 = open('cluster-Al.xyz','a+')#these two output file is a collection of xyz format data, which is centered by Al and B atoms, respectively
fout2 = open('cluster-B.xyz','a+')
for i in range(1,na+1):
    nbl2xyz(neighborlist(i,1,n,fin),fout1,fout2)
fout1.close()
fout2.close()
