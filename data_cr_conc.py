# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 02:52:06 2017

@author:
"""

import os
import pickle

GENES = map(lambda s: s.replace('.seq', ''), os.listdir("genes"))
genedict = {}

for i, gene in enumerate(GENES):
    genedict[i] = gene

mutated_genes = []

for i in range(len(genedict)):
    s = "enc_gene_" + genedict[i] + ".csv"
    print (s)
    if os.path.isfile(s):
        mutated_genes.append(i)
        
def sampling(g_dict, m_genes):
    non_mutated = {}
    concat = {}
    for i in xrange(len(g_dict)):
        f = open("genes/" + str(g_dict[i]) + ".csv", "rb")
        non_mutated[i] = f.read()
        
    for i in mutated_genes:

        current_gene = g_dict[i]
        f = pickle.load(open(s = "enc_gene_" + g_dict[i] + ".csv"),"rb")
        
        new_first_part = str()
        new_second_part = str()
        
        for x in xrange(g_dict[i]):
            new_first_part += str(x)
            
        for x in range(g_dict[i] + 1, len(g_dict)):
            new_second_part += str(x)
                
        n = f.shape[0]

        for j in xrange(n):
            new_seq = new_first_part
            
            if ((f[0,j] == ( str(g_dict[i]) + "_sub" ) ) ):  
                
                new_seq += str(f[1,j])
                new_seq += new_second_part
            
            concat[g_dict[i]].append(new_seq)
        
    pickle.dump(concat, open("concat" + ".csv", "wb"))
    
