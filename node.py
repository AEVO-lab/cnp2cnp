import random
import sys
import argparse
import os
from subprocess import Popen, PIPE
from shutil import copyfile
import math


class Node:

	id = None
	data = None
	children = []
	parent = None

	def __init__(self):
		pass

	def get_copy(self):
		n = Node()
		n.id = self.id
		n.data = self.data

		for c in children:
			copy = c.get_copy()
			copy.parent = n
			n.children.append(copy)

		return n

	def add_child(self):
		n = Node()
		n.parent = self
		children.append(n)
		return n

	def is_leaf(self):
		return (len(self.children) == 0)

	def is_root(self):
		return self.parent == None

	@staticmethod
	def get_random_binary_tree(nb_leaves):
		tree = Node.get_random_binary_tree_rec(nb_leaves)

		nodes = tree.get_postorder_nodes()

		i = 1
		for n in nodes:
			if n.is_leaf():
				n.id = i
				i += 1

		return tree

	def to_newick(self, print_data = False):

		if print_data:
			data_str = str(self.data).replace("[", "").replace("]", "").replace(",", "-").replace(" ", "")

		if self.is_leaf():
			if print_data:
				strret = data_str
				return str(self.id) + "__" + strret
			else:
				return str(self.id)
		else:

			chstr = []
			for c in self.children:
				chstr.append(c.to_newick(print_data))


			strret = "(" + ','.join(chstr) + ")"

			if print_data:
				strret += data_str

			return strret

	@staticmethod
	def get_random_binary_tree_rec(nb_leaves):
		if nb_leaves == 1:
			return Node()

		nleft = random.randint(1,nb_leaves - 1)

		leftch = Node.get_random_binary_tree_rec(nleft)
		rightch = Node.get_random_binary_tree_rec(nb_leaves - nleft)

		n = Node()
		n.children = [leftch, rightch]
		leftch.parent = n
		rightch.parent = n

		return n

	def get_postorder_nodes(self, list = None):
		if list == None:
			list = []

		for c in self.children:
			c.get_postorder_nodes(list)

		list.append(self)

		return list
