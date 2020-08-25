# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import sys

# create a new 'Legacy VTK Reader'
a020614_1_11_144_skeletonvtk = LegacyVTKReader(FileNames=['/data/proj/mitochondrialNetworks/Viana2020MendeleyDataset/020614_1_11_144_skeleton.vtk'])

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
a020614_1_11_000_skeletonvtk = GetActiveSource()

# save data
SaveData('/data/proj/mitochondrialNetworks/Viana2020MendeleyDataset/test/testfile1.csv', proxy=a020614_1_11_000_skeletonvtk, PointDataArrays=['Intensity', 'Width'])

print("finished")

sys.exit()
