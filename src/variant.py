#!/usr/bin/env python

class Variant:
    def __init__(self, name, chrom, pos, orig, alts, probs):
        self.name = name#7
        self.chrom = chrom#0
        self.pos = pos#1
        self.orig = orig#2
        self.alts = alts#3
        self.probs = probs#4
        self.num_alts = len(alts)

    def add_alt(self, alt, prob):
        self.alts.append(alt)
        self.probs.append(prob)
        self.num_alts += 1
