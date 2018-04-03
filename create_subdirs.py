#!/usr/bin/env python

'''
Creating subdirs for assembly
--------------------------------------------------------------------------------------

BinMate - Metagenome binning pipeline
J .Frank - j.frank@science.ru.nl
P .Boersma - p.boersma@science.ru.nl

Description
Quick script for seperating fastq files needed for assembly.

Usage
create_subdirs.py -i <path/to/dir/>

--------------------------------------------------------------------------------------
'''
import argparse

def sub_dir(args):
	print (args.input)
	dir_list = []
	for file in args.input:
		file1 = file.split("/")[1]
		print (file)
	print (dir_list)

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Assembly Statistics Script')
        parser.add_argument('-i', '--input',
                            type=str,
                            nargs="+",
                            required=True,
                            help='Directory to fastq files')
        args = parser.parse_args()

        sub_dir(args)

