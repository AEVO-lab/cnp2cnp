import random
import sys
import argparse
import os
from subprocess import Popen, PIPE
from shutil import copyfile
import math
import numpy

import node
from cnpsolver import CNPSolver


class GenomeSimulator:

	prob_dup = 0.5
	max_duplength_ratio = 0.5
	max_losslength_ratio = 0.5
	max_duplength = 500
	max_losslength = 500

	prob_stop = 0.05
	prob_removelast = 0.1

	def __init__(self):
		pass

	def get_default_genome(self, nbgenes):
		return list(range(nbgenes))



	'''
	source: a sequence of arbitrary elements, each int representing a character
	'''
	def simulate_events(self, source, nb_events, print_them_all = False):

		dest = source
		for i in range(nb_events):
			dest = self.simulate_event(dest)
			if print_them_all:
				print(dest)

		return dest

	def simulate_event(self, source_p):

		if len(source_p) == 0:
			return source_p


		source = list(source_p)  #copy because we modify

		#coin determines dup or loss
		coin = random.random()

		

		is_dup = False
		if float(coin) <= float(self.prob_dup):
			is_dup = True


		'''
		max_len = self.max_duplength
		max_len_r = self.max_duplength_ratio
		if not is_dup:
			max_len_r = self.max_losslength_ratio
			max_len = self.max_losslength

		lenratio = random.random() * max_len_r


		#k = 1 + int(lenratio * (len(source) - 1))	#k is the length
		#k = min(max_len, k)
		k = 1
		prob_stop = 0.1
		done_k = False
		while not done_k:
			coin = random.random()
			
			if float(coin) > float(prob_stop)  and k < len(source):
				adv = 1 #int(random.random() * 20)
				k += adv
			else:
				done_k = True
		#print("k=" + str(k))

		s = int(random.random() * (len(source) - k))
		t = (min(len(source), s + k))

		#print("k = " + str(k) + " s = " + str(s) + " t = " + str(t))
		'''
		
		prob_stop = self.prob_stop

		if is_dup:
			s = int(random.random() * (len(source) - 1))
			t = s + 1
			done = False
		

			while not done:
				if t == len(source) or t - s >= self.max_duplength:
					done = True
				else:
					coin = random.random()
					if (float(coin) <= float(prob_stop)):
						done = True
					else:
						t += 1
			copy = source[s : t]



			#dups are tandem
			dest_tmp = source[0:s] + copy + source[s:]
			
			#check that nothing is above 9
			maxgid = 0
			for v in dest_tmp:
				if v > maxgid:
					maxgid = v
			cnp = CNPSolver.get_cnp_from_genome(dest_tmp, maxgid + 1)
			isok = True

			#THIS BELOW WAS TO ACCOMODATE MEDICC, BUT NO, NOT ANYMORE
			#for c in cnp:
			#	if c > 4:
			#		isok = False

			if isok:
				dest = dest_tmp
			else:
				dest = source

		else:
			s = int(random.random() * (len(source) - 1))
			done = False
			prob_removelast = self.prob_removelast
			delcount = 0

			while not done:
				if s == len(source) or delcount >= self.max_losslength:
					done = True
				else:
					maxgid = 0
					for v in source:
						if v > maxgid:
							maxgid = v
					cnp = CNPSolver.get_cnp_from_genome(source, maxgid + 1)
					coin = random.random()
					if (float(coin) <= float(prob_stop)):
						done = True
					else:
						cnt = cnp[source[s]]
						if cnt == 1:	#if we are about to remove the last guy
							coin = random.random()
							#print(str(coin) + "   vs   " + str(prob_removelast))
							if float(coin) <= float(prob_removelast):
								del source[s]
								delcount += 1
							else:
								done = True
						else:
							del source[s]
							delcount += 1
				

			dest = source


		return dest


	def simulate_tree_evolution(self, tree, root_genome, min_events, max_events):
		tree.data = root_genome

		for c in tree.children:
			nbevents = min_events + int(random.random() * (max_events - min_events))
			r = numpy.random.normal(1, 1)
			nbevents = int(max(1, round(nbevents * r)))
			


			dest = self.simulate_events(root_genome, nbevents)
			self.simulate_tree_evolution(c, dest, min_events, max_events)
