# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:00:56 2017

@author:
"""

import os
import pandas as pd
import re
from collections import Counter

GENES = map(lambda s: s.replace('.seq', ''), os.listdir("genes"))
print(GENES)


del_re = re.compile(r'c\.(?P<pos>\d+)del(?P<nucl>[ATGC])')  # regexp for deletion
ins_re = re.compile(r'c\.(?P<pos>\d+)ins(?P<nucl>[ATGC])')  # regexp for insertion
sub_re = re.compile(r'c\.(?P<pos>\d+)(?P<old_nucl>[ATGC])>(?P<new_nucl>[ATGC])')  # regexp for substitution
mutations_regexs = (del_re, ins_re, sub_re)

def apply_mutation(seq, mutation):
    new_seq = None
    msg = lambda pos, n_exp, n_real, mut: "Nucleotides not matching in " + mut + \
            "! Expected " + n_exp + " on pos " + str(pos) + ", but actual is " + n_real
    
    m = re.search(sub_re, mutation)
    if m:
        pos = int(m.group('pos')) - 1  # because in Python arrays begin with zero
        old_n = m.group('old_nucl')
        new_n = m.group('new_nucl')
        if seq[pos] != old_n:
            raise Exception(msg(pos, old_n, seq[pos], mutation))
        return seq[:pos] + new_n + seq[pos+1:], "sub"
    
    m = re.search(del_re, mutation)
    if m:
        pos = int(m.group('pos')) - 1  # because in Python arrays begin with zero
        n = m.group('nucl')
        if seq[pos] != n:
            raise Exception(msg(pos, n, seq[pos], mutation))
        return seq[:pos] + seq[pos+1:], "del"
        
    m = re.search(ins_re, mutation)
    if m:
        pos = int(m.group('pos')) - 1  # because in Python arrays begin with zero
        n = m.group('nucl')
        if seq[pos] != n:
            raise Exception(msg(pos, n, seq[pos], mutation))
        return seq[:pos] + n + seq[pos:], "ins"
        
    return None, None

def gen_seq_arr(GENE):
    with open("genes/" + GENE + ".seq") as f_in:
        gene_seq = ''.join(map(str.rstrip, f_in.readlines()))
        
    csv_f = pd.read_csv("genes_mutations/" + GENE + "_mutations.csv", delimiter='\t')
    csv_f = csv_f[csv_f['Primary site'] == 'lung']
    if len(gene_seq) != csv_f['Gene CDS length'].unique()[0]:
        print ("Gene " + GENE + " has length " + str(len(gene_seq)) + \
            ", but mutations are expected to have " + str(csv_f['Gene CDS length'].unique()[0]))
        return None, None
    seq_arr = []
    c = Counter()
    for mutation in csv_f['Mutation CDS']:
        new_seq, m_type = apply_mutation(gene_seq, mutation)
        if new_seq:
            seq_arr.append(str(new_seq))
            c[m_type] += 1
    return seq_arr, c

for GENE in GENES:
    seq_arr, c = gen_seq_arr(GENE)
    if seq_arr:
        print("Gene " + GENE + " has " + str(len(seq_arr)) + " mutant variants: {}".format(c))
        s = "enc_gene_" + GENE + ".csv"
        gene_out = open(s, "w")
        gene_out.write(str(seq_arr))
