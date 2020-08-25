# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:37:22 2020

@author: ekoslover
"""
import glob
from paraview.simple import *

# get list of files to process
dir = '/data/proj/mitochondrialNetworks/Viana2020MendeleyDataset/';

filenames = glob.glob(dir+'*skeleton.vtk')

for filename in filenames:
    print(filename)
    
    outfile = filename.replace(".vtk",".txt")
    

    # load structure into paraview
    # create a new 'Legacy VTK Reader'
    skeletonvtk = LegacyVTKReader(FileNames=[filename])
    skeletonvtk = GetActiveSource()
    
    SaveData(outfile, proxy=skeletonvtk, PointDataArrays=[])