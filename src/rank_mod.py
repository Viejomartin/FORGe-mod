#!/usr/bin/env python

'''
Rank a set of variants for inclusion in a graph genome, from highest to lowest priority
'''

import sys
import argparse
import iohelp
from util import *
import pandas as pd
import numpy as np
from pysam import FastaFile
import resource
import os
import time

AUTHOR = 'Martín Prado'
VERSION = '0.0.1 (modified)'

class VarRanker1:
    def __init__(self, genome, variants, r, phasing, max_v):
        self.genome = genome
        self.chrom_lens = dict()
        for chrom, seq in genome.items():
            self.chrom_lens[chrom] = len(seq)

        self.variants = variants
        self.num_v = variants.shape[0]
        self.r = r

        self.phasing = phasing
        self.hap_parser = iohelp.HaplotypeParser(phasing) if phasing else None

        self.max_v_in_window = max_v

        self.h_ref = None
        self.h_added = None

        self.wgt_ref = None
        self.wgt_added = None

        self.curr_vars = None
        self.freqs = {}

    def rank(self, method, out_file):
        ordered = None
        ordered_blowup = None
        print(method)
        if method == 'popcov':
            ordered = self.rank_pop_cov()
        elif method == 'popcov-blowup':
            ordered = self.rank_pop_cov(True)
        elif method == 'hybrid':
            ordered, ordered_blowup = self.rank_hybrid()

        if ordered is not None:
            df = pd.DataFrame(columns = ['a', 'b'])
            df['a'] = self.variants.loc[ordered].chrom
            df['b'] = self.variants.loc[ordered].pos +1
            df.to_csv(out_file,header=False,index=False,line_terminator = '\t')
            # with open(out_file, 'w') as f:
            #     f.write('\t'.join([str(self.variants.chrom[i]) + ',' for i in ordered]))
        if ordered_blowup:
            df = pd.DataFrame(columns = ['a', 'b'])
            df['a'] = self.variants.loc[ordered_blowup].chrom
            df['b'] = self.variants.loc[ordered_blowup].pos +1
            df.to_csv(out_file,header=False,index=False,line_terminator = '\t')
            # with open(out_file+'.blowup', 'w') as f:
            #     f.write('\t'.join([self.variants.chrom[i] + ',' + str(self.variants.pos[i]+1) for i in ordered_blowup]))

    def rank_pop_cov(self, with_blowup=False, threshold=1/3):
        '''
        Rank variants using the hybrid ranking method.
        
        with_blowup: If true, add blowup penalty for neighboring variants
        threshold: Blowup penalty for neighboring variants, between 0 and 1. A smaller value is a stricter penalty
        '''
        
        if with_blowup:
            varsdf = self.variants
            var_wgts = (self.variants.iloc[:,[4,8,9]].sum(axis = 1)).sort_index()
            posdf = varsdf.pos
            ffirst = np.vectorize(
                lambda i: 
                    (posdf <= posdf[i]-self.r).idxmin()
                    )

            flast = np.vectorize(
                lambda i: 
                    (posdf >= posdf[i]+self.r).idxmax()
                    )
            first = np.maximum(0,np.fromfunction(ffirst,(self.num_v,), dtype=int))
            last = np.minimum(self.num_v,np.fromfunction(flast,(self.num_v,), dtype=int))-1
            last[last == -1] = self.num_v-1
            
            neighbors = last-first
            penalty = np.repeat(threshold,self.num_v)
            var_wgts = pd.Series(var_wgts* (penalty**neighbors).copy()).sort_values(ascending=False)
            ordered = var_wgts.index
        else:
            # Variant weight is the sum of frequencies of alternate alleles
            var_wgts = (self.variants.iloc[:,[4,8,9]].sum(axis = 1)).sort_index()
            var_wgts = var_wgts.sort_values(ascending=False)
            ordered = var_wgts.index

        return ordered


def go(args):
    start_time = time.time()
    x1 = 'Memory usage start: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if args.window_size:
        r = args.window_size
    else:
        r = 10
    if args.output:
        o = args.output
    else:
        o = 'ordered.txt'
    if args.prune:
        max_v = args.prune
    else:
        max_v = r
    #genome = iohelp.read_genome(args.reference, args.chrom)
    genome1 = iohelp.read_genome1(args.reference, args.chrom)
    x2 = 'Memory usage after  genome created: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    vars1 = iohelp.parse_1ksnp1(args.vars)
    x3 = 'Memory usage after  variants created: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ranker = VarRanker1(genome1, vars1, r, args.phasing, max_v)
    ranked = ranker.rank(args.method, o)
    x4 = 'Memory usage after  ranking created: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    x5 = "Running time: %s seconds" % (time.time() - start_time)
    oram = 'ram-'+args.method+'-mod.txt'
    with open(oram, 'w') as f:
        f.writelines([x1,'\n',x2,'\n',x3,'\n',x4,'\n',x5])

    


if __name__ == '__main__':

    if '--version' in sys.argv:
        print('ERG v' + VERSION)
        sys.exit(0)

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--method', type=str, required=True,
        help='Variant ranking method. Currently supported ranking methods: [popcov | popcov-blowup | hybrid]\n\'hybrid\' will produce hybrid ranking files both with and without blowup avoidance,')
    parser.add_argument('--reference', type=str, required=True, 
        help='Path to fasta file containing reference genome')
    parser.add_argument("--vars", type=str, required=True,
        help="Path to 1ksnp file containing variant information")
    parser.add_argument('--chrom', type=str, required=True,
        help="Name of chromosome from reference genome to process. If not present, process all chromosomes.")
    parser.add_argument('--window-size', type=int,
        help="Radius of window (i.e. max read length) to use. Larger values will take longer. Default: 10")
    parser.add_argument('--pseudocontigs', action="store_true", help='Rank pseudocontigs rather than SNPs')
    parser.add_argument('--phasing', type=str, required=False,
        help="Path to file containing phasing information for each individual")
    parser.add_argument('--output', type=str, required=False,
        help="Path to file to write output ranking to. Default: 'ordered.txt'")
    parser.add_argument('--prune', type=int, required=False,
        help='In each window, prune haplotypes by only processing up to this many variants. We recommend including this argument when ranking with the hybrid strategy for window sizes over 10.')

    args = parser.parse_args(sys.argv[1:])
    go(args)
