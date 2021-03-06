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

def apply_mutation(seq, mutation, path_status=0):
    new_seq = None
    msg = lambda pos, n_exp, n_real, mut: "Nucleotides not matching in " + mut + \
            "! Expected " + n_exp + " on pos " + str(pos) + ", but actual is " + n_real
    
    if path_status == 0:
        suffix = "_?"
    elif path_status == 1:
        suffix = "_P"
    elif path_status == 2:
        suffix = "_N"

    m = re.search(sub_re, mutation)
    if m:
        pos = int(m.group('pos')) - 1  # because in Python arrays begin with zero
        old_n = m.group('old_nucl')
        new_n = m.group('new_nucl')
        if seq[pos] != old_n:
            raise Exception(msg(pos, old_n, seq[pos], mutation))
        return seq[:pos] + new_n + seq[pos+1:], "sub" + suffix
    
    m = re.search(del_re, mutation)
    if m:
        pos = int(m.group('pos')) - 1  # because in Python arrays begin with zero
        n = m.group('nucl')
        if seq[pos] != n:
            raise Exception(msg(pos, n, seq[pos], mutation))
        return seq[:pos] + seq[pos+1:], "del" + suffix
        
    m = re.search(ins_re, mutation)
    if m:
        pos = int(m.group('pos')) - 1  # because in Python arrays begin with zero
        n = m.group('nucl')
        if seq[pos] != n:
            raise Exception(msg(pos, n, seq[pos], mutation))
        return seq[:pos] + n + seq[pos:], "ins" + suffix
        
    return None, None

def gen_seq_arr(GENE):
    with open("genes/" + GENE + ".seq") as f_in:
        gene_seq = ''.join(map(str.rstrip, f_in.readlines()))
    seq_arr = [(GENE, gene_seq)]
        
    csv_f = pd.read_csv("genes_mutations/" + GENE + "_mutations.csv", delimiter='\t')
    csv_f = csv_f[csv_f['Primary site'] == 'lung']
    if len(gene_seq) != csv_f['Gene CDS length'].unique()[0]:
        print ("Gene " + GENE + " has length " + str(len(gene_seq)) + \
            ", but mutations are expected to have " + str(csv_f['Gene CDS length'].unique()[0]))
        return None, None
    c = Counter()
    for i, row in csv_f.iterrows():
        #print i
        #print csv_f['FATHMM prediction']
        mutation = row['Mutation CDS']
        if row['FATHMM prediction'] == 'PATHOGENIC':
            path_status = 1
        elif row['FATHMM prediction'] == 'NEUTRAL':
            path_status = 2
        else:
            path_status = 0
        new_seq, m_type = apply_mutation(gene_seq, mutation, path_status)
        if new_seq:
            seq_arr.append((GENE + "_" + m_type, str(new_seq)))
            c[m_type] += 1
    return seq_arr, c

for GENE in GENES:
    seq_arr, c = gen_seq_arr(GENE)
    if seq_arr:
        print("Gene " + GENE + " has " + str(len(seq_arr)) + " mutant variants: {}".format(c))
        s = "enc_gene_" + GENE + ".csv"
        with open(s, "w") as gene_out:
            for seq_tup in seq_arr:
                m_type, seq = seq_tup
                gene_out.write(m_type + '\t' + seq + '\n')
