import random
import sys
import argparse
import os
from subprocess import Popen, PIPE
from shutil import copyfile
from genomesimulator import GenomeSimulator
import math



indir = "results_rcg"

for fname in os.listdir(indir):
	if fname.endswith(".csv"):
		ff = open(indir + "/" + fname)
		nb = 0
		ttl_us = 0
		ttl_us_flat = 0
		ttl_euc = 0
		ttl_zzs = 0
		for line in ff:
			line = line.replace("\n", "")
			if line != "":
				nb += 1
				if nb > 1:
					pz = line.split(",")
					us = float(pz[6])
					us_flat = float(pz[7])
					euc = float(pz[8])
					zzs = float(pz[9])
					ttl_us += us
					ttl_us_flat += us_flat
					ttl_euc += euc
					ttl_zzs += zzs
		
		avg_us = ttl_us/float(nb)
		avg_us_flat = ttl_us_flat/float(nb)
		avg_euc = ttl_euc/float(nb)
		avg_zzs = ttl_zzs/float(nb)
		print(fname + "," + str(avg_us) + "," + str(avg_us_flat) + "," + str(avg_zzs) + "," + str(avg_euc))
			
			
