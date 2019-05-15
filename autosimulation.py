import random
import sys
import argparse
import os
from subprocess import Popen, PIPE
from shutil import copyfile
from genomesimulator import GenomeSimulator
import math


#note that this file is called autosimulation.py and not autostimulation.py
outdir = "results_dbl"
runs = 100


'''
gs = GenomeSimulator()
source = gs.get_default_genome(10)


gs.prob_dup = 0.5
dest = gs.simulate_events(source, 100)

print(source)
print(dest)

sys.exit()
'''


	



#isolate number of genes
for g in [10,50,100,500]:
	outx = outdir + "/sim_genes_" + str(g)
	if not os.path.exists(outx):
	    os.mkdir(outx)

	cmd = "python cnp_simulator.py -m simulate_trees -n 100 -g " + str(g) + " -d " + outx + " -r " + str(runs) + " -p 0.5 -e 5:10 -o " + outdir + "/genes_" + str(g) + ".csv"
	print("EXEC " + cmd)
	os.system(cmd)
	#sys.exit()


#isolate number of leaves
for l in [10,50,100,500]:
	outx = outdir + "/sim_leaves_" + str(l)
	if not os.path.exists(outx):
	    os.mkdir(outx)
	cmd = "python cnp_simulator.py -m simulate_trees -n " + str(l) + " -g 100 -d " + outx + " -r " + str(runs) + " -p 0.6 -e 5:10 -o " + outdir + "/leaves_" + str(l) + ".csv"
	print("EXEC " + cmd)
	os.system(cmd)





#isolate duprate
for d in [0.25, 0.5, 0.75]:
	outx = outdir + "/sim_dr_" + str(d)
	if not os.path.exists(outx):
	    os.mkdir(outx)
	cmd = "python cnp_simulator.py -m simulate_trees -n 100 -g 100 -d " + outx + " -r " + str(runs) + " -p " + str(d) + " -e 5:10 -o " + outdir + "/duprate_" + str(d) + ".csv"
	print("EXEC " + cmd)
	os.system(cmd)

#isolate nbevents
for e in ["2:4", "5:10", "20:40"]:
	outx = outdir + "/sim_ev_" + e.replace(":", "_")
	if not os.path.exists(outx):
	    os.mkdir(outx)
	cmd = "python cnp_simulator.py -m simulate_trees -n 100 -g 100 -d " + outx + " -r " + str(runs) + " -p 0.6 -e " + e + " -o " + outdir + "/events_" + e.replace(":", "_") + ".csv"
	print("EXEC " + cmd)
	os.system(cmd)
