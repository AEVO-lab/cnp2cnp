import argparse
from ete2 import Tree
import sys

parser = argparse.ArgumentParser()
parser.add_argument('treefile1')
parser.add_argument('treefile2')

args = parser.parse_args()

t1 = Tree(args.treefile1)
t2 = Tree(args.treefile2)

rf, rf_max, x1, x2, x3, x4, x5 = t1.robinson_foulds(t2, unrooted_trees=True)
print(str(rf) + "/" + str(rf_max))
