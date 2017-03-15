# Quick script to easily extract the names of the tip labels from a given tree

import argparse
import traceback
import warnings
from Bio import Phylo


def main():
	try:
		parser = argparse.ArgumentParser(description='Retrieve the names of tip labels in a phylogenetic tree')
		parser.add_argument(
			'-i',
			'--infile',
			action='store',
			help='The treefile to extract names from.')
		parser.add_argument(
			'-f',
			'--format',
			action='store',
			default='newick',
			help='The file format of the treefile. Default expects NEWICK.')
		args = parser.parse_args()

	except:
		print('An exception occurred when parsing the arguments. Check your provided options, or see --help.')



	tree = Phylo.read(args.infile, args.format)
	for leaf in tree.get_terminals():
		print leaf.name

if __name__ == "__main__":
	main()
