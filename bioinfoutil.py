import random
import sys
import argparse
import os
import os.path
from subprocess import Popen, PIPE
from shutil import copyfile
import subprocess
import math


class BioinfoUtil:

	@staticmethod
	def write_matrix_to_phylip(outfile, names, matrix):
		strout = str(len(matrix)) + "\n"

		for i in range(len(matrix)):
			strout += names[i].ljust(10) + " "
			for j in range(len(matrix[i])):
				if j != 0:
					strout += " "
				strout += str(matrix[i][j])
			strout += "\n"
		f = open(outfile, 'w')
		f.write(strout)
		f.close()

	@staticmethod
	def get_rf_dist(treefile1, treefile2, format = None):
		cmd = "python rfdist.py " + treefile1 + " " + treefile2
		#if format != None:
		#	cmd += " " + str(format)
		proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
		(out, err) = proc.communicate()
		
		pz = out.split("/")
		return [int(pz[0]), int(pz[1])]


	@staticmethod
	def run_phylip_nj(distfile, outfile, workdir = ''):

				
		if os.path.exists(os.path.join(workdir, "infile")):
  			os.remove(os.path.join(workdir, "infile"))
		if os.path.exists(os.path.join(workdir, "outfile")):
  			os.remove(os.path.join(workdir, "outfile"))
		if os.path.exists(os.path.join(workdir, "outtree")):
  			os.remove(os.path.join(workdir, "outtree"))

		copyfile(distfile, os.path.join(workdir, "infile"))

		oldcwd = os.getcwd()
		
		os.chdir(workdir)
		

		pr = Popen(['phylip', 'neighbor'], stdout=PIPE, stderr=PIPE, stdin=PIPE)

		pr.stdin.write("y\n")

		pr.wait()

		os.chdir(oldcwd)
		copyfile(os.path.join(workdir, "outtree"), outfile)

		#pr.stdin.write(distfile + "\n")

        ## say yes


		#cmd = "phylip neighbor -datafile " + distfile + " -outfile " + outfile
		#os.system(cmd)
