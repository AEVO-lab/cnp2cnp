import random
import sys
import argparse
import os
from os.path import join
from subprocess import Popen, PIPE
from shutil import copyfile
from genomesimulator import GenomeSimulator
import math


#note that this file is called autosimulation.py and not autostimulation.py

runs = 100

default_nbleaves = 100
default_nbgenes = 100
default_probdup = 0.5
default_events = "5:10"
doerror = True

for prst in [0.01, 0.05, 0.1]:
	for prls in [0.25, 0.5, 0.75, 1]:
		default_probstop = prst
		default_probloss = prls
		error_rate = 1


		outdir = "results3_s" + str(default_probstop).replace(".", "") + "l" + str(default_probloss).replace(".", "") + "x" + str(error_rate).replace(".", "")

		#medicc_param = " -q $HOME/projects/medicc/medicc/"




		if not os.path.exists(outdir):
		    os.mkdir(outdir)

		#-n number of leaves
		#-g number of genes (positions)
		#-p prob dup
		#-e events per branch




		#isolate number of genes
		for g in [100,250,10,50]: 
			outx = outdir + "/sim_genes_" + str(g)
			if not os.path.exists(outx):
			    os.mkdir(outx)

			cmd = "python cnp_simulator.py -m simulate_trees -n " + str(default_nbleaves) + " -g " + str(g) + " -d " + outx + " -r " + str(runs) + " -p " + str(default_probdup) + " -e " + default_events + " -o " + outdir + "/genes_" + str(g) + ".csv"
			cmd += " -s " + str(default_probstop) + " -l " + str(default_probloss)
			cmd += " -x " + str(error_rate)
			if not doerror:
				cmd += " -z 1"
			#cmd += medicc_param
			print("EXEC " + cmd)
			os.system(cmd)




		#isolate number of leaves
		for l in [10, 50, 100]:
			outx = outdir + "/sim_leaves_" + str(l)
			if not os.path.exists(outx):
			    os.mkdir(outx)
			cmd = "python cnp_simulator.py -m simulate_trees -n " + str(l) + " -g " + str(default_nbgenes) + " -d " + outx + " -r " + str(runs) + " -p " + str(default_probdup) + " -e " + default_events + " -o " + outdir + "/leaves_" + str(l) + ".csv"
			cmd += " -s " + str(default_probstop) + " -l " + str(default_probloss)
			cmd += " -x " + str(error_rate)
			#cmd += medicc_param
			if not doerror:
				cmd += " -z 1"

			print("EXEC " + cmd)
			os.system(cmd)





		#isolate duprate
		for d in [0.25, 0.5, 0.75]:
			outx = outdir + "/sim_dr_" + str(d)
			if not os.path.exists(outx):
			    os.mkdir(outx)
			cmd = "python cnp_simulator.py -m simulate_trees -n " + str(default_nbleaves) + " -g " + str(default_nbgenes) + " -d " + outx + " -r " + str(runs) + " -p " + str(default_probdup) + " -e " + default_events + " -o " + outdir + "/duprate_" + str(d) + ".csv"
			cmd += " -s " + str(default_probstop) + " -l " + str(default_probloss)
			cmd += " -x " + str(error_rate)

			if not doerror:
				cmd += " -z 1"


			print("EXEC " + cmd)
			#cmd += medicc_param
			os.system(cmd)

		#isolate nbevents
		for e in ["2:4", "5:10", "20:40"]:
			outx = outdir + "/sim_ev_" + e.replace(":", "_")
			if not os.path.exists(outx):
			    os.mkdir(outx)
			cmd = "python cnp_simulator.py -m simulate_trees -n " + str(default_nbleaves) + " -g " + str(default_nbgenes) + " -d " + outx + " -r " + str(runs) + " -p " + str(default_probdup) + " -e " + e + " -o " + outdir + "/events_" + e.replace(":", "_") + ".csv"
			cmd += " -s " + str(default_probstop) + " -l " + str(default_probloss)
			cmd += " -x " + str(error_rate)


			if not doerror:
				cmd += " -z 1"

			print("EXEC " + cmd)
			#cmd += medicc_param
			os.system(cmd)
