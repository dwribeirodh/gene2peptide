#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 17:40:02 2020

@author: danielribeiro
"""

import numpy as np
import statistics
import random
import time

#Let A = 0
#Let C = 1
#Let G = 2
#Let T = 3


def nuc_seq(nsites):
    init = tuple([0.0,3.0,2.0])
    sequence = [init]
    for i in range(nsites):
        l = []
        for a in range(3):
            n = random.randint(0,4)
            l = np.append(l, n)
        l = tuple(l)
        sequence.append(l)
    return sequence
        

def gene2peptide(sequence):
    #pdb.set_trace()
    codon_table = {tuple([3,3,3]):"F", tuple([3,3,1]):"F", tuple([3,3,0]):"L", tuple([3,3,2]):"L",
                   tuple([3,1,3]):"S", tuple([3,1,1]):"S", tuple([3,1,0]):"S", tuple([3,1,2]):"S",
                   tuple([3,0,3]):"Y", tuple([3,0,1]):"Y", tuple([3,0,0]):"_", tuple([3,0,2]):"_",
                   tuple([3,2,3]):"C", tuple([3,2,1]):"C", tuple([3,2,0]):"_", tuple([3,2,2]):"W",
                   tuple([1,3,3]):"L", tuple([1,3,1]):"L", tuple([1,3,0]):"L", tuple([1,3,2]):"L",
                   tuple([1,1,3]):"P", tuple([1,1,1]):"P", tuple([1,1,0]):"P", tuple([1,1,2]):"P",
                   tuple([1,0,3]):"H", tuple([1,0,1]):"H", tuple([1,0,0]):"Q", tuple([1,0,2]):"Q",
                   tuple([1,2,3]):"R", tuple([1,2,1]):"R", tuple([1,2,0]):"R", tuple([1,2,2]):"R",
                   tuple([0,3,3]):"I", tuple([0,3,1]):"I", tuple([0,3,0]):"I", tuple([0,3,2]):"M",
                   tuple([0,1,3]):"T", tuple([0,1,1]):"T", tuple([0,1,0]):"T", tuple([0,1,2]):"T",
                   tuple([0,0,3]):"N", tuple([0,0,1]):"N", tuple([0,0,0]):"K", tuple([0,0,2]):"K",
                   tuple([0,2,3]):"S", tuple([0,2,1]):"S", tuple([0,2,0]):"R", tuple([0,2,2]):"R",
                   tuple([2,3,3]):"V", tuple([2,3,1]):"T", tuple([2,3,0]):"T", tuple([2,3,2]):"T",
                   tuple([2,1,3]):"A", tuple([2,1,1]):"A", tuple([2,1,0]):"A", tuple([2,1,2]):"A",
                   tuple([2,0,3]):"D", tuple([2,0,1]):"D", tuple([2,0,0]):"E", tuple([2,0,2]):"E",
                   tuple([2,2,3]):"G", tuple([2,2,1]):"G", tuple([2,2,0]):"G", tuple([2,2,2]):"G"}
    
    peptide = []
    for codon in sequence:
        if codon in codon_table:
            if codon == (3,0,0) or codon == (3,2,0) or codon == (3,0,2):
                break
            else:
                peptide += codon_table[codon]
    return peptide

def stats_sequence(nsites, nseq):
    
    nseqs = []
    print("Generating gene sequences...")
    for i in range(nseq):
        rand_seq = nuc_seq(nsites)
        protein = gene2peptide(rand_seq)
        nseqs.append(protein)
    
    print("Calculating AA chain length...")
    length = np.array([])
    for sequence in nseqs:
        sequence_length = len(sequence)
        length = np.append(length, sequence_length)
    
    print("Getting stats...")
    mean = np.mean(length)
    mode = statistics.mode(length)
    median = np.median(length)
    
    return mean, mode, median
start_time = time.time()
if __name__ == "__main__":
    nsites = 1500
    nseqs = 60000
    seq_mean, seq_mode, seq_median = stats_sequence(nsites, nseqs)
    print("--- %s seconds ---" % (time.time() - start_time))