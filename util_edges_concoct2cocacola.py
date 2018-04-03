#!/usr/bin/env python

"""
Convert CONCOCT lfile to COCACOLA edge_list
--------------------------------------------------------------------------------

BinMate - Metagenome binning pipeline
J .Frank - j.frank@science.ru.nl

Usage
util_edges_concoct2cocacola.py -i <input_linkage_concoct.tsv> -o <output_edges_cocacola.tsv>

Parses CONCOCT bam_to_linkage.py output (linkage between contigs)
to COCACOLA format: <contig_i> <contig_n> <links>
For each contig pair, this script will select either
inward, outward or inline paired-ed connections depending on which has highest number of links

--------------------------------------------------------------------------------
"""

import argparse
import csv

parser = argparse.ArgumentParser(description="util_edges_concocnt2cocacola")

parser.add_argument("-i", "--input", type=str, required=True,
                            help="Input: CONCOCT linkage file")
parser.add_argument("-o", "--output", type=str, required=True,
                            help="Output COCACOLA edges file")

args = parser.parse_args()


with open(args.input, "r") as input:
  links = csv.reader(input, delimiter="\t")
  size = len(next(links))
          
  # Select linkage types
  inward = list(range(2, size, 6))
  outward = list(range(3, size, 6))
  inline = list(range(4, size, 6))

  out_edges = open(args.output, "w")
  edges = csv.writer(out_edges, delimiter="\t")

  # For each contig pair, write most link type with highest number of links to file (inward, outward or inline)
  for dat in links:
    edges.writerow([dat[0], dat[1], max(sum([int(dat[i]) for i in inward]), sum([int(dat[i]) for i in outward]), sum([int(dat[i]) for i in inline]))])

out_edges.close()
