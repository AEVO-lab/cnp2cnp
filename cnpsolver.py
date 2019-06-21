import random
import sys
import argparse
import os
from subprocess import Popen, PIPE
from shutil import copyfile
import math


class CNPSolver:
	@staticmethod
	def get_nb_flat_intervals(u, v):

		count = 0
		prev_w = 0
		for i in range(len(u)):
			delta = u[i] - v[i]

			if delta != prev_w:
				count += 1
			prev_w = delta

		if prev_w == 0 and count > 0:
			count -= 1

		return count


	'''
	static method
	genes are numbered from 1 through nbgenes
	'''
	@staticmethod
	def get_cnp_from_genome(genome, nbgenes):
		cnp = [0] * nbgenes

		for i in range(len(genome)):
			cnp[genome[i]] += 1

		return cnp

	'''
	TODO: use something clever like numpy to do this
	'''
	@staticmethod
	def get_cnp_difference(u, v):
		w = []
		for i in range(len(u)):
			w.append(u[i] - v[i])
		return w

	@staticmethod
	def get_comparable_cnps(u, v):
		up = []
		vp = []
		cnt_zerodiff = 0

		up.append(1)
		vp.append(1)	#does not change the distance, but prevents extreme zeros

		for i in range(len(u)):
			if u[i] > 0:
				up.append(u[i])
				vp.append(v[i])
			elif u[i] == 0 and v[i] != 0:
				cnt_zerodiff += 1

		up.append(1)
		vp.append(1)	#does not change the distance, but prevents extreme zeros
		return [up, vp, cnt_zerodiff]


	@staticmethod
	def get_ZZS_distance(src, dest):
		N = 0
		for i in range(len(src)):
			if src[i] > N:
				N = src[i]
			if dest[i] > N:
				N = dest[i]

		M = []	#dimensions are number of positions x N
		for i in range(len(src)):
			if dest[i] != 0:
				subarr = [0] * (N + 1)
				
				for d in range(N + 1):
					if d >= src[i]:
						subarr[d] = 99999
					elif d < src[i] - dest[i]:
						subarr[d] = 99999
					elif i == 0:
						subarr[d] = d - (src[i] - dest[i])
					else:
						#main recurrence
						
						#prev position is closest non-zero to the left
						prev = 1
						for k in range(i):
							if dest[k] > 0:
								prev = k
						curmin = 99999
						for dp in range(N + 1):
							aid = d - (src[i] - dest[i])
							aprev = dp - (src[prev] - dest[prev])
							sum = M[prev][dp] + max(d - dp, 0) + max(aid - aprev, 0)
							
							q = 0
							for j in range(prev + 1, i):
								if dest[j] == 0 and src[j] > q:
									q = src[j]
							sum += max(q - max(d, dp), 0)
							if sum < curmin:
								curmin = sum
						subarr[d] = curmin
				M.append(subarr)
			else:
				M.append([])
	
		mincost = 99999
		for d in range(N + 1):
			if M[len(src) - 1][d] < mincost:
				mincost = M[len(src) - 1][d]

		return mincost

	@staticmethod
	def read_fasta_cnp_file(filename, single_digit_mode = False):
		lst = []
		cnps = []
		f = open(filename)
		for line in f:
			line = line.replace("\n", "")
			if line.startswith(">"):
				lst.append(line.replace(">", ""))
			elif line != "":
				sz = []
				if single_digit_mode:
					for c in range(len(line)):
						sz.append(line[c])			
				else:
					sz = line.split(",")
				cnp = []
				for val in sz:
					cnp.append(int(val))
				cnps.append(cnp)

		return [lst, cnps]

	@staticmethod
	def get_distance_matrix(cnps, flats = False, use_dbl = False, zzs = False):
		#init
		the_matrix = []
		for i in range(len(cnps)):
            		the_matrix.append([0] * len(cnps))

		for i in range(len(cnps)):
			for j in range(i + 1, len(cnps)):
				up = cnps[i]
				vp = cnps[j]
				tmp = CNPSolver.get_comparable_cnps(up, vp)
				u = tmp[0]
				v = tmp[1]
				cntz = tmp[2]

				if zzs:
				    score = CNPSolver.get_ZZS_distance(u, v)
				    the_matrix[i][j] = score		
				    the_matrix[j][i] = score
				elif not flats:
				    evs = CNPSolver.get_approximate_events(u, v, use_dbl)
				    the_matrix[i][j] = len(evs) #+ cntz
				    the_matrix[j][i] = len(evs) #+ cntz
				else:
				    nbflat = CNPSolver.get_nb_flat_intervals(u, v)
				    the_matrix[i][j] = nbflat
				    the_matrix[j][i] = nbflat
		return the_matrix

	@staticmethod
	def get_euclidean_distance_matrix(cnps):
		#init
		matrix = []
		for i in range(len(cnps)):
			matrix.append([0] * len(cnps))

		for i in range(len(cnps)):
			for j in range(i + 1, len(cnps)):
				u = cnps[i]
				v = cnps[j]

				sum = 0
				for k in range(len(u)):
					sum += (u[k] - v[k])*(u[k] - v[k])
				dist = math.sqrt(float(sum))

				matrix[i][j] = dist
				matrix[j][i] = dist

		return matrix

	@staticmethod
	def get_approximate_events(u_param, v, use_dbl = False):

		u = list(u_param)	#make a copy, otherwise u gets modified

		done = False
		events = []

		while not done:

			ev = CNPSolver.get_best_event(u, v)


			if ev == None:
				done = True
			else:
				events.append(ev)
				s = ev[0]
				t = ev[1]
				b = ev[2]

				maxb = b
				if b > 0:
					for p in range(s, t + 1):
						maxb = u[p]
				if use_dbl:
					ev[2] = min(b, maxb)
					b = min(b, maxb)


				strout = str(u) + "(" + str(ev) + ")="

				for i in range(s, t + 1):
					if u[i] != 0:
						u[i] = max(u[i] + b, 0)

				strout += str(u)


		return events

	@staticmethod
	def get_first_nonzero(w):
		for i in range(len(w)):
			if w[i] != 0:
				return i
		return None

	@staticmethod
	def get_last_nonzero(w):
		for i in reversed(range(len(w))):
			if w[i] != 0:
				return i
		return None

	@staticmethod
	def get_best_event(u, v):

		w = CNPSolver.get_cnp_difference(u, v)

		s = CNPSolver.get_first_nonzero(w)
		if s == None:
			return None
		t = CNPSolver.get_last_nonzero(w)




		if s == t:
			return [s, t, -w[s]]

		one_flat = True
		for i in range(s + 1, t + 1):
			if w[i] != w[i-1]:
				one_flat = False

		if one_flat:
			return [s, t, -w[s]]



		#ok, so all trivial cases are taken care of...


		pos_per_delta = {}



		#this is the naive brute-force implementation
		for i in range(s, t + 1):

			if i == s:
				pos_per_delta[w[i]] = s
			else:
				if i == t:
					delta_next = w[i]
				else:
					delta_next = w[i] - w[i+1]

				if delta_next != 0 and u[i] != 0:

					if delta_next in pos_per_delta:
						pos = pos_per_delta[delta_next]

						can_do_it = True
						for j in range(pos, i + 1):
							if u[j] <= delta_next and v[j] > 0:
								can_do_it = False

						if can_do_it:
							return [pos, i, -delta_next]

				delta_prev = w[i] - w[i - 1]

				if delta_prev != 0:
					pos_per_delta[delta_prev] = i

		#if we didn't find anything, then just reduce the first non-zero flat interval
		start_interval = s
		for i in range(s + 1, t + 1):
			if w[i] != w[i-1]:
				return [s, i - 1, -w[s]]

		#if we make it here, all of w is flat, and was verified above
		print("ERROR 12: how did you get here?")
		sys.exit()
