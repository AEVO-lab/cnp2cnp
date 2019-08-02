import random
import sys
import argparse
import os
from os.path import join
from subprocess import Popen, PIPE
from shutil import copyfile
from shutil import rmtree
import math
import ntpath

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
parser.add_argument('-w', "--workdir", default='work', help='Work dir for temp files in mode infer_tree_medicc')
parser.add_argument('-r', "--nbruns", default='100', help='number of runs (of trees to simulate)')
parser.add_argument('-e', "--nbevents", default='10:20', help='min:max of eventson a tree branch (format is a:b)')
parser.add_argument('-p', "--probdup", default='0.75', help='probability of duplication (beware, <.7 tends to remove every gene eventually)')
parser.add_argument('-x', "--errorrate", default='0', help='.  0 = no error is introduced.')

parser.add_argument('-q', "--mediccdir", default='', help='directory that contains medicc.py')

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
			
if args.mode == 'fasta_to_medicc':
	infile = args.infile
	[names, cnps] = CNPSolver.read_fasta_cnp_file(infile)
	
	strout = ""
	
	#medicc needs a diploid
	strout += ">diploid\n"
	for i in cnps[0]:
		strout += "1"
	strout += "\n"

	for i in range(len(names)):
		strout += ">" + names[i] + "\n"
		vals = cnps[i]
		for v in vals:
			if int(v) >= 10:
				print("ERROR: a copy-number of sequence " + names[i] + " is greater than 9.  Medicc will not run.")
				sys.exit()
			strout += str(v)
		strout += "\n"
	if args.outfile == '':
		print(strout)
	else:
		f = open(args.outfile, 'w')
		f.write(strout)
		f.close()


if args.mode == 'infer_tree_medicc':
	infile = args.infile
	
	if not os.path.exists(args.workdir):
	    os.mkdir(args.workdir)

	workdir = join(args.workdir, "medicc_tmp")

	if not os.path.exists(workdir):
	    os.mkdir(workdir)

	#clean the workdir	
	for the_file in os.listdir(workdir):
		file_path = os.path.join(workdir, the_file)

		if os.path.isfile(file_path):
			os.unlink(file_path)

	cmd = "python cnp_simulator.py -m fasta_to_medicc -i " + infile + " -o "
	

	majorfile = join(workdir, ntpath.basename(infile) + ".major_chr1.fasta")
	minorfile = join(workdir, ntpath.basename(infile) + ".minor_chr1.fasta")

	

	#medicc handles and requires a minor and major allele file.  For our tests, we make them equal
	#print("EXEC " + cmd + majorfile)	
	os.system(cmd + majorfile)
	#print("EXEC " + cmd + minorfile)
	os.system(cmd + minorfile)


	if not os.path.isfile(majorfile) :
		print("*** FAILED TO GENERATE MAJOR FILE ***")
		sys.exit()
	#print("*** SUCCEEDED IN GENERATING MAJOR FILE ***")

	#medicc also requires a desc file
	descfilename = join(workdir, ntpath.basename(infile) + ".desc.txt")
	fdesc = open(descfilename, 'w')
	fdesc.write("chrom1 " + ntpath.basename(majorfile) + " " + ntpath.basename(majorfile))
	fdesc.close()

	#run medicc	-s prevents ancestor reconstructions, -n says 'don't resolve alleles'
	mediccout_dir = join(workdir, "medicc.out")
	if os.path.isdir(mediccout_dir) :
		rmtree(mediccout_dir)
	cmd = "python " + join(args.mediccdir, "medicc.py") + " " + descfilename + " " + mediccout_dir + " -s"
	print("EXEC " + cmd)
	os.system(cmd)

	medicc_treefile = join(workdir, "medicc.out", "tree_fitch_nc.new")
	copyfile(medicc_treefile, args.outfile)
	

if args.mode == 'infer_tree' or args.mode == 'infer_tree_euclidean' or args.mode == 'infer_tree_flat' or args.mode == 'infer_tree_zzs':
	infile = args.infile
	#print(infile)
	[names, cnps] = CNPSolver.read_fasta_cnp_file(infile)

	if args.mode == 'infer_tree':
		matrix = CNPSolver.get_distance_matrix(cnps)
	elif args.mode == 'infer_tree_flat':
		matrix = CNPSolver.get_distance_matrix(cnps, flats = True)
	elif args.mode == 'infer_tree_zzs':
		matrix = CNPSolver.get_distance_matrix(cnps, zzs = True)
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

	line = "TREE,RF_HEUR,RF_FLAT,RF_EUCLID,RF_ZZS,MAXRF,NORMRF_HEUR,NORMRF_FLAT,NORMRF_EUCLID,NORMRF_ZZS"
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


		#we don't support medicc anymore, it doesn't work well
		#test with medicc, or at least try
		#inferredtreefilename_medicc = outdir + "/tree" + str(r) + ".inferred_medicc.newick"

		#if os.path.isfile(inferredtreefilename_medicc):
		#	os.unlink(inferredtreefilename_medicc)
		#cmd = "python cnp_simulator.py -m infer_tree_medicc -i " + fastafilename + " -w " + join(outdir, "work") + " -q " + args.mediccdir + " -o " + inferredtreefilename_medicc
		#print("EXEC " + cmd)
		#os.system(cmd)

		#if not os.path.isfile(inferredtreefilename_medicc):
		#	print("*** MEDICC MADE AN ERROR, SKIPPING ***")
		#	continue

		#[rfdist_medicc, maxrf_medicc] = BioinfoUtil.get_rf_dist(treefilename, inferredtreefilename_medicc, format=1)


		
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


		#infer with zzs distance
		inferredtreefilename_zzs = outdir + "/tree" + str(r) + ".inferred_zzs.newick"
		cmd = "python cnp_simulator.py -m infer_tree_zzs -i " + fastafilename + " -o " + inferredtreefilename_zzs
		#print("EXEC " + cmd)
		os.system(cmd)

		[rfdist_zzs, maxrf_zzs] = BioinfoUtil.get_rf_dist(treefilename, inferredtreefilename_zzs)
		#print("RF=" + str(rfdist) + "/" + str(maxrf))

		line = treefilename + "," + str(rfdist_us) + "," + str(rfdist_flat) + "," + str(rfdist_them) + "," + str(rfdist_zzs) + "," + str(maxrf_us)
		line += "," + str(float(rfdist_us)/float(maxrf_us)) + "," + str(float(rfdist_flat)/float(maxrf_us)) + "," + str(float(rfdist_them)/float(maxrf_them)) + "," + str(float(rfdist_zzs)/float(maxrf_zzs))

		print(line)
		if args.outfile != "":
			fout.write(line + "\n")

		os.system("rm -rf ~/.local/share/Trash/*")

	if args.outfile != "":
		fout.close()



