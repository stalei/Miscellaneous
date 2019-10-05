#  © Shahram Talei @ 2019 The University of Alabama
# yt is required:
# https://github.com/yt-project

import yt
import numpy as np
from yt.analysis_modules.halo_finding.api import *
from yt.analysis_modules.halo_analysis.api import *
from os import environ
environ['CFLAGS'] = "-I"+np.get_include()

import pyximport; pyximport.install()
#import particle_ops
import argparse


import tempfile
import shutil
import os
import sys


unit_base = {'length': (1.0, 'Mpc'),
    'velocity': (1.0, 'km/s'),
    'mass': (1.0, 'Msun')}

def _IsContaminated(halo,ds, Rv):#snap):
	#halo
	dds = halo.halo_catalog.data_ds
	center = dds.arr([halo.quantities["particle_position_%s" % axis] \
	for axis in "xyz"])
	my_id = halo.quantities['particle_identifier']
	print("Halo %d" % (my_id))
	#particles
	count=0
	boundary=np.zeros((2,3))
	ad = ds.all_data()
	coordinatesDM = ad[("Halo","Coordinates")]
	coordinatesLowRes = ad[("Bndry","Coordinates")]
	DMCount=len(coordinatesDM)
	LowResCount=len(coordinatesLowRes)
	print("number of High resolution particles:%d"%DMCount)
	print("number of Low resolution particles:%d"%LowResCount)
	for i in range(0,3):
		boundary[0,i]=np.min(coordinatesDM[:,i]) # (min max, x y z)
		boundary[1,i]=np.max(coordinatesDM[:,i])
	print("high resolution particles are in box:")
	print(boundary)
	#for i in range(0,LowResCount): # this for is super slow!! 
	#	x=coordinatesLowRes[i,0]
	#	y=coordinatesLowRes[i,1]
	#	z=coordinatesLowRes[i,2]
	#	if(x>boundary[0,0] and x<boundary[1,0] and y>boundary[0,1] and y<boundary[1,1] and z>boundary[0,2] and z<boundary[1,2]):
	#		print("contamination at: %f,%f,%f"%x%y%z)
	#		count+=1
	xLR=coordinatesLowRes[:,0]
	yLR=coordinatesLowRes[:,1]
	zLR=coordinatesLowRes[:,2]
	
	xLR2=xLR[(xLR>boundary[0,0]) & (xLR<boundary[1,0])]
	yLR2=yLR[(xLR>boundary[0,0]) & (xLR<boundary[1,0])]
	zLR2=zLR[(xLR>boundary[0,0]) & (xLR<boundary[1,0])]
	
	xLR3=xLR2[(yLR2>boundary[0,1]) & (yLR2<boundary[1,1])]
	yLR3=yLR2[(yLR2>boundary[0,1]) & (yLR2<boundary[1,1])]
	zLR3=zLR2[(yLR2>boundary[0,1]) & (yLR2<boundary[1,1])]
	
	
	xLR4=xLR3[(zLR3>boundary[0,2]) & (zLR3<boundary[1,2])]
	yLR4=yLR3[(zLR3>boundary[0,2]) & (zLR3<boundary[1,2])]
	zLR4=zLR3[(zLR3>boundary[0,2]) & (zLR3<boundary[1,2])]
	
	r2=((xLR4-center[0])**2.+(yLR4-center[1])**2.+(zLR4-center[2])**2.)
	r=np.sqrt(r2)
	# should be Rv instead of 0.17 but I have to check their unit and convert them.
	rContamination=r[r<0.17]
	count=len(rContamination)
	print("Halo %d at:" % (my_id))
	print(center)
	print("is contaminated with: %d low resolution particles at δr="%count)
	print(rContamination)
	#return count

# Hot to run?
# python contamination.py snap_264
# snapshot is in gadget binary format. 
# Snapshot is zoom-in file produced by music initial condition. 
#High resolution area is type "Dark matter" and low resolution is "Boundary particle"
#With few changes it should work for other types of the simulations.

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("snap", type=str)
	
	args = parser.parse_args()
	ds = yt.load(args.snap)
	if ds is None:
		print ("Error, couldn't load your parameter files!")
		sys.exit(1)
	
	hc = HaloCatalog(data_ds=ds, finder_method='hop', finder_kwargs={ "dm_only": False , "ptype":"Halo"} )
	hc.create(save_halos=True)
	ad = hc.halos_ds.all_data()
	masses = ad['particle_mass'][:]
	Rvir = ad['virial_radius'][:]
	add_callback("IsContaminated", _IsContaminated)
	hc.add_callback("IsContaminated", ds, Rvir)
	hc.create(save_halos=True)
