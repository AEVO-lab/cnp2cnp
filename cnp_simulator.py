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

parser = argparse.ArgumentParser(description='Get the number of deletion and duplication events from a source CNP to a target CNP.')
parser.add_argument('-m', "--mode", default='', help='mode of execution')
parser.add_argument('-d', "--outdir", default='.', help='output directory')
parser.add_argument('-i', "--infile", default='', help='input file')
parser.add_argument('-o', "--outfile", default='', help='output file')
parser.add_argument('-n', "--nbleaves", default='10', help='number of leaves')
parser.add_argument('-g', "--nbgenes", default='10', help='number of genes')
parser.add_argument('-r', "--nbruns", default='100', help='number of runs (of trees to simulate)')
parser.add_argument('-e', "--nbevents", default='10:20', help='min:max of eventson a tree branch (format is a:b)')
parser.add_argument('-p', "--probdup", default='0.75', help='probability of duplication (beware, <.7 tends to remove every gene eventually)')

args = parser.parse_args()




if args.mode == 'concat_medicc_majors':
	indir = args.infile
	seqs = {}
	for file in os.listdir(indir):
		if file.endswith("major.fasta"):
			fullfile = os.path.join(indir, file)
			[names, cnps] = CNPSolver.read_fasta_cnp_file(fullfile, single_digit_mode = True)
			for i in range(len(names)):
				if names[i] not in seqs:
					seqs[names[i]] = []
				seqs[names[i]] += cnps[i]
	for s in seqs:
		print(">" + s)
		print(seqs[s])
			


if args.mode == 'infer_tree' or args.mode == 'infer_tree_euclidean' or args.mode == 'infer_tree_flat':
	infile = args.infile
	#print(infile)
	[names, cnps] = CNPSolver.read_fasta_cnp_file(infile)

	if args.mode == 'infer_tree':
		matrix = CNPSolver.get_distance_matrix(cnps)
	elif args.mode == 'infer_tree_flat':
		matrix = CNPSolver.get_distance_matrix(cnps, flats = True)
	else:
		matrix = CNPSolver.get_euclidean_distance_matrix(cnps)


	distfile = args.outfile + ".dist"
	BioinfoUtil.write_matrix_to_phylip(distfile, names, matrix)
	BioinfoUtil.run_phylip_nj(distfile, args.outfile)



if args.mode == 'simulate_trees':
	nbruns = int(args.nbruns)
	outdir = args.outdir

	strout = ""

	if args.outfile != "":
		fout = open(args.outfile, 'w')

	line = "TREE,RF_HEUR,RF_FLAT,RF_EUCLID,MAXRF,NORMRF_HEUR,NORMRF_FLAT,NORMRF_EUCLID"
	print(line)
	if args.outfile != "":
		fout.write(line + "\n")


	for r in range(nbruns):

		nb_genes = int(args.nbgenes)
		nb_leaves = int(args.nbleaves)
		gs = GenomeSimulator()
		source = gs.get_default_genome(nb_genes)

		#we do a few initial dups, as otherwise the genome tends to lose everything
		gs.prob_dup = 1
		source = gs.simulate_events(source, 10)
		tree = Node.get_random_binary_tree(nb_leaves)

		minev = int(args.nbevents.split(":")[0])
		maxev = int(args.nbevents.split(":")[1])


		gs.prob_dup = args.probdup
		
		gs.simulate_tree_evolution(tree, source, minev, maxev)

		#write CNPs in a modified fasta format: characters are counts and are comma separated
		strfa = ""
		for n in tree.get_postorder_nodes():
			#n.data = CNPSolver.get_cnp_from_genome(n.data, nb_genes)
			if n.is_leaf():
				n.id = "taxon_" + str(n.id)
				strfa += ">" + n.id + "\n"

				cnp = CNPSolver.get_cnp_from_genome(n.data, nb_genes)
				for i in range(len(cnp)):
					if i != 0:
						strfa += ","
					strfa += str(cnp[i])
				strfa += "\n"

		fastafilename = outdir + "/cnps" + str(r) + ".fa"
		fastaf = open(fastafilename, 'w')
		fastaf.write(strfa)
		fastaf.close()

		treefilename = outdir + "/tree" + str(r) + ".newick"
		treefile = open(treefilename,'w')
		treefile.write(tree.to_newick(False) + ";")
		treefile.close()

		treefilename2 = outdir + "/tree" + str(r) + ".detailed.newick"
		treefile2 = open(treefilename2,'w')
		treefile2.write(tree.to_newick(True) + ";")
		treefile2.close()


		#test with our heuristic approx
		inferredtreefilename = outdir + "/tree" + str(r) + ".inferred.newick"
		cmd = "python cnp_simulator.py -m infer_tree -i " + fastafilename + " -o " + inferredtreefilename
		#print("EXEC " + cmd)
		os.system(cmd)

		[rfdist_us, maxrf_us] = BioinfoUtil.get_rf_dist(treefilename, inferredtreefilename)

		#test with the naive flat count
		inferredtreefilename_flat = outdir + "/tree" + str(r) + ".inferred_flat.newick"
		cmd = "python cnp_simulator.py -m infer_tree_flat -i " + fastafilename + " -o " + inferredtreefilename_flat
		#print("EXEC " + cmd)
		os.system(cmd)

		[rfdist_flat, maxrf_flat] = BioinfoUtil.get_rf_dist(treefilename, inferredtreefilename_flat)

		#infer with euclidean distance
		inferredtreefilename_euc = outdir + "/tree" + str(r) + ".inferred_euc.newick"
		cmd = "python cnp_simulator.py -m infer_tree_euclidean -i " + fastafilename + " -o " + inferredtreefilename_euc
		#print("EXEC " + cmd)
		os.system(cmd)

		[rfdist_them, maxrf_them] = BioinfoUtil.get_rf_dist(treefilename, inferredtreefilename_euc)
		#print("RF=" + str(rfdist) + "/" + str(maxrf))

		line = treefilename + "," + str(rfdist_us) + "," + str(rfdist_flat) + "," + str(rfdist_them) + "," + str(maxrf_us)
		line += "," + str(float(rfdist_us)/float(maxrf_us)) + "," + str(float(rfdist_flat)/float(maxrf_us)) + "," + str(float(rfdist_them)/float(maxrf_them))

		print(line)
		if args.outfile != "":
			fout.write(line + "\n")

	if args.outfile != "":
		fout.close()



