#!/usr/bin/env python

'''
Assembly statistics report
--------------------------------------------------------------------------------

BinMate - Metagenome binning pipeline
J .Frank - j.frank@science.ru.nl

Description
Quick & dirty script to assess characteristics of an assembly.
Based on a script (assemblathon.pl) by Keith Bradnam.

Usage
qc_asm_stats.py -i <contigs.fasta> -o <out_dir> --report

Input
- assembled contigs/scaffolds (fasta)
Output
- statistics report (.txt)
- data table (tsv)

Dependencies
- Biopython
- Python package: Pandas

Citations

Cock, P. J. et al. 
Biopython: freely available Python tools for computational molecular biology and bioinformatics. 
Bioinformatics 25, 1422-1423, doi:10.1093/bioinformatics/btp163 (2009).

--------------------------------------------------------------------------------
'''

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC
import os
import sys
import datetime
import time
import pandas
import argparse


parser = argparse.ArgumentParser(description='Assembly Statistics Script')

parser.add_argument('-i', '--file', type=str, required=True, help='Contigs/scaffolds FASTA file')
parser.add_argument('-l', '--length', type=str, required=False, default=int(0), help='Filter contigs/scaffolds on length (bp)')
parser.add_argument('-n', '--number', type=int, required=False, help='Return top n largest contigs/scaffolds (int)')
parser.add_argument('-r', '--report', action='store_true', help='Create assembly statistics report')
parser.add_argument('-d', '--datatable', action='store_true', help='Write datatable')
parser.add_argument('-f', '--fasta', action='store_true', help='Write FASTA file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output directory')

args = parser.parse_args()
startTimer = time.time()


def anaylse_fasta(fasta, lenfilt):
    '''
    Generate contig/scaffold data
    '''

    data = []
    with open(fasta) as input:
        for id, seq in SimpleFastaParser(input):
            if len(seq) >= int(lenfilt):
                record = []
                seq = seq.upper()

                record.append(id)
                record.append(len(seq))
                record.append(GC(seq))
                record.append(seq.count('A'))
                record.append(seq.count('T'))
                record.append(seq.count('C'))
                record.append(seq.count('G'))
                record.append(seq.count('N'))
                data.append(record)

    header = ['contig_id', 'size', 'gc', 'a', 't', 'c', 'g', 'n_ambiguous']
    df = pandas.DataFrame(data, columns=header).sort_values('size', ascending=False)
    
    if args.number is not None:
        if df.shape[0] >= args.number:
            df = df.head(args.number)
        else:
            args.number = df.shape[0]
            
    return(df)


def report(df):
    '''
    Calculate statistics and generate report
    '''

    tot_scf     = df.shape[0]
    tot_bp      = df['size'].sum()
    max_size    = df['size'].max()
    min_size    = df['size'].min()
    mean_scf    = df['size'].mean()
    median_scf  = df['size'].median()
    kb1         = (df['size'] > 1000).sum()
    kb1_p       = kb1*100.0/tot_scf
    kb10        = (df['size'] > 10000).sum()
    kb10_p      = kb10*100.0/tot_scf
    kb100       = (df['size'] > 100000).sum()
    kb100_p     = kb100*100.0/tot_scf
    mb1         = (df['size'] > 1000000).sum()
    mb1_p       = mb1*100.0/tot_scf
    mean_scf_gc = df['gc'].mean()
    a_bp        = (df['a'].sum())*100.0/tot_bp
    t_bp        = (df['t'].sum())*100.0/tot_bp
    c_bp        = (df['c'].sum())*100.0/tot_bp
    g_bp        = (df['g'].sum())*100.0/tot_bp
    n           = df['n_ambiguous'].sum()
    n_bp        = n*100.0/tot_bp

    # Calculate N50
    l50 = 0
    sum = 0
    for size in df.sort_values('size', ascending=False)['size']:
        if sum >= tot_bp*0.5:
            n50 = df.iloc[l50]['size']
            break
        else:
            sum += size
            l50 += 1

    # Print report to file
    date = datetime.datetime.now().strftime('%d-%m-%Y - %H:%M %Z')
    if int(args.length) != 0:
        output = '{}/asm-stats_{}_gt{}.txt'.format(args.output+"asm_stats", filename, args.length)
    elif args.number is not None:
        output = '{}/asm-stats_{}_top{}.txt'.format(args.output+"asm_stats", filename, args.number)
    else:
        output = '{}/asm-stats_{}.txt'.format(args.output+"asm_stats", filename)
    report = open(output, 'w')

    original_stdout = sys.stdout
    sys.stdout = report

    print('Assembly Statistics Script\n')
    print('BinMate - Metagenome binning pipeline\n')
    print('J .Frank - j.frank@science.ru.nl\n\n')
    print(date)
    print('File: {}'.format(args.file))
    
    if int(args.length) != 0:
       print('Filter: excluded contigs < {} bp'.format(args.length))
    if args.number is not None:
        print('Filter: stats for top {} largest contigs/scaffolds'.format(args.number))

    print('\nNumber of scaffolds:       {:,}'.format(tot_scf))
    print('Total assembled bp:        {:,}'.format(tot_bp))
    print('Longest scaffold:          {:,}'.format(max_size))
    print('Shortest scaffold:         {:,}'.format(min_size))
    print('Mean scaffold size:        {:,}'.format(round(mean_scf,2)))
    print('Median scaffold size:      {:,}'.format(round(median_scf,2)))
    print('No. scaffolds > 1 Kb:      {:,}\t({}%)'.format(kb1, round(kb1_p,2)))
    print('No. scaffolds > 10 Kb:     {:,}\t({}%)'.format(kb10, round(kb10_p,2)))
    print('No. scaffolds > 100 Kb:    {:,}\t({}%)'.format(kb100, round(kb100_p,2)))
    print('No. scaffolds > 1 Mb:      {:,}\t({}%)'.format(mb1, round(mb1_p,2)))
    print('N50 scaffold length:       {:,}'.format(n50))
    print('L50 scaffold count:        {:,}'.format(l50))
    print('Mean scaffold GC:          {}%'.format(round(mean_scf_gc,2)))
    print('Nuc A:                     {}%'.format(round(a_bp,2)))
    print('Nuc T:                     {}%'.format(round(t_bp,2)))
    print('Nuc C:                     {}%'.format(round(c_bp,2)))
    print('Nuc G:                     {}%'.format(round(g_bp,2)))
    print('Nuc ambiguous (Ns):        {:,}\t({}%)'.format(n, round(n_bp,2)))

    print('\nProcessing time: {} seconds'.format(round(time.time() - startTimer, 2)))

    sys.stdout = original_stdout
    report.close()


def save_data(df):
    if int(args.length) != 0:
        output  = '{}/ctg-stats_{}_gt{}.tsv'.format(args.output+"/ctg_stats", filename, args.length)
    elif args.number is not None:
        output = '{}/ctg-stats_{}_top{}.tsv'.format(args.output+"/ctg_stats", filename, args.number)
    else:
        output = '{}/ctg-stats_{}.tsv'.format(args.output+"/ctg_stats", filename)

    df.to_csv(output, sep='\t', float_format='%.2f', quoting=3, index=False)

def save_fasta(df):
    if int(args.length) != 0:
        output  = '{}/{}_gt{}.fasta'.format(args.output, filename, args.length)
    elif args.number is not None:
        output = '{}/{}_top{}.fasta'.format(args.output, filename, args.number)

    record_dict = SeqIO.index(args.file, 'fasta')
    output = open(output, 'w')
    [SeqIO.write(record_dict[seq.strip('>\n').split()[0]], output, 'fasta') for seq in df['contig_id']]

if not any ((args.report, args.datatable)):
  print('No options selected. Exit')
  sys.exit(0)
else:
  filename = os.path.splitext(os.path.basename(args.file))[0]
  data = anaylse_fasta(args.file, args.length)

if args.report:
    report(data)
if args.datatable:
    save_data(data)
if args.fasta:
    save_fasta(data)
