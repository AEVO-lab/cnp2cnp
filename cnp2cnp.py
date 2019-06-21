import random
import sys
import argparse
import os
from subprocess import Popen, PIPE
from shutil import copyfile
import math

from node import Node
from cnpsolver import CNPSolver
from genomesimulator import GenomeSimulator
from bioinfoutil import BioinfoUtil

#write 
#> python cnp2cnp.py --help to see available options
parser = argparse.ArgumentParser(description='Get the number of deletion and duplication events from a source CNP to a target CNP.')
parser.add_argument('-m', "--mode", default='dist', help='Mode of execution.  Put "dist" to compute the distance between two CNPs and output it to stdout, "matrix" to get a distance matrix from CNPs in a fasta file and output it in a phylip file.')
parser.add_argument('-i', "--infile", default='', help='Input file.  A fasta file for modes "dist" (with two CNPs) or "matrix".  CNP values must be comma separated.')
parser.add_argument('-o', "--outfile", default='', help='Output file, or stdout if none specified.  Only valid for mode "matrix".')
parser.add_argument('-d', "--distance", default='any', help='The distance to use.  Possible values are "any": use improved approximation from Cordonnier and Lafond, "flat": count the number of flat intervals, "dbl": use approximation for "any", but forbids counts to be above their double in a single event, "euclidean" for the euclidean distance, "zzs" to use the Zeira-Zehavi-Shamir distance.')

args = parser.parse_args()





if args.mode == 'dist':
	'''
	mode "dist"
	Compute the distance between two CNPs.  Expected input file is in fasta format and has 2 CNPs.
	Example of fasta file: 
	> my_first_cnp
	1,2,0,3,2,1
	> my_second_cnp
	1,2,3,0,1,1

	The positions with a 0 in the first CNP will be removed.
	The computed distance will be written to stdout.
	'''
	if args.infile == "":
		print("You must specify an input file.  It should be a fasta file with two CNPs.")
		sys.exit()

	infile = args.infile
	[names, cnps] = CNPSolver.read_fasta_cnp_file(infile)
	
	if len(names) > 2:
		print("Your fasta file has more than 2 sequences.  Did you want '-m matrix'?")
		sys.exit()
	if len(names) < 2:
		print("Your fasta file has less than 2 sequences.  Why?")
		sys.exit()
	
	tmp = CNPSolver.get_comparable_cnps(cnps[0], cnps[1])
	u = tmp[0]
	v = tmp[1]


	if args.distance == "any" or args.distance == "dbl":
		p_use_dbl = (args.distance == "dbl")
		evs = CNPSolver.get_approximate_events(u, v, use_dbl = p_use_dbl)
		print(len(evs))
	elif args.distance == "flat":
		nbflat = CNPSolver.get_nb_flat_intervals(u, v)
		print(nbflat)
	elif args.distance == "euclidean":
		matrix = CNPSolver.get_euclidean_distance_matrix(cnps)
		print(matrix[0][1])
	elif args.distance == "zzs":
		d = CNPSolver.get_ZZS_distance(u, v)
		print(d)

elif args.mode == 'matrix':
	'''
	mode "matrix"
	Compute a distance matrix between all pairs of sequences in a given fasta file (see above for format).
	In case of "any", "dbl" or "flat" distances, the distance is not symmetric, but the output matrx is symmetric.  
	The program will compute d(u, v) if u is before v in the fasta file, or d(v, u) otherwise.
	The output is written to a phylip file, destined to be used by the "phylip neighbor" program.
	'''


	if args.infile == "":
		print("You must specify an input file.  It should be a fasta file with two or more CNPs.")
		sys.exit()

	if args.outfile == "":
		print("You must specify an output file for the distance matrix (in phylip format)")
		sys.exit()

	infile = args.infile
	[names, cnps] = CNPSolver.read_fasta_cnp_file(infile)
	
	if args.distance == "any" or args.distance == "dbl":
		p_use_dbl = (args.distance == "dbl")
		matrix = CNPSolver.get_distance_matrix(cnps, flats = False, use_dbl = p_use_dbl)
	elif args.distance == "flat":
		matrix = CNPSolver.get_distance_matrix(cnps, flats = True)
	elif args.distance == "zzs":
		matrix = CNPSolver.get_distance_matrix(cnps, flats = True, zzs = True)
	elif args.distance == "euclidean":
		matrix = CNPSolver.get_euclidean_distance_matrix(cnps)

	BioinfoUtil.write_matrix_to_phylip(args.outfile, names, matrix)
		





		


