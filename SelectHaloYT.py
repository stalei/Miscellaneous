
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




if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("snap", type=str)
	args = parser.parse_args()
	ds = yt.load(args.snap)
	if ds is None:
		print ("Error, couldn't load your parameter files!")
		sys.exit(1)
	#halos = HaloFinder(ds)
	#ind = halos[len(halos)-1]["particle_index"] # list of particles IDs in this halo
	#mass = halos[0]["particle_mass"]/#"particle_mass"
	#print(halos)
	#mass_solar=mass/5.02785e-34
	#print(mass_solar)
	# Load the rockstar data files
	#halos = yt.load('halos_0.0.bin')
	# Create temporary directory for storing files
	#tmpdir = tempfile.mkdtemp()
	# Instantiate a catalog using those two parameter files
	#hc = HaloCatalog(data_ds=ds, halos_ds=halos, output_dir=os.path.join(tmpdir, 'halo_catalog'))
	hc = HaloCatalog(data_ds=ds, finder_method='hop')
	# Filter out less massive halos
	hc.add_filter("quantity_value", "particle_mass", ">", 1.2e12, "Msun")
	hc.add_filter("quantity_value", "particle_mass", "<", 1.4e12, "Msun")
	# Save the profiles
	hc.add_callback("save_profiles", storage="virial_profiles", output_dir="profiles")
	hc.create()