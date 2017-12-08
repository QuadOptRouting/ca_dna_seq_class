# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 02:52:06 2017

@author:
"""

import os
import pickle
import pandas as pd
import csv

GENES = map(lambda s: s.replace('.seq', ''), os.listdir("genes"))
genedict = {}

for i, gene in enumerate(GENES):
    genedict[i] = gene

mutated_genes = {}
ii = 0

for i in range(len(genedict)):
    s = "enc_gene_" + genedict[i] + ".csv"
    print (s)
    if os.path.isfile(s):
        mutated_genes[ii] = genedict[i]
        ii += 1


print mutated_genes

ii = 0


def sampling(g_dict):
    non_mutated = {}
    concat = {}
    gene_tables = {}
    f = open("conc.csv", 'wb')
    for i in xrange(len(g_dict)):
        gene_tables[i] = pd.read_csv("enc_gene_" + str(g_dict[i]) + ".csv", sep='\t', header=None, names=['name', 'sequence'])
        non_mutated[i] = gene_tables[i].iloc[0]['sequence']

    # for i in xrange(len(g_dict)):
    #     temp = gene_tables[i]['sequence']
    #     tmp = 0.
    #     for s in temp:
    #         tmp += len(s)
    #     print tmp/temp.shape[0]

    for i in mutated_genes:
        current_gene = g_dict[i]
        print current_gene
        concat[current_gene] = []
        
        new_first_part = str()
        new_second_part = str()
        
        for x in xrange(i):
            new_first_part += str(non_mutated[x])
        for x in range(i + 1, len(g_dict)):
            new_second_part += str(non_mutated[x])
        n = gene_tables[i].shape[0]
        flg = True
        for j in xrange(1, n):
            new_seq = new_first_part
            if flg:
                flg = False
            if gene_tables[i].iloc[j]['name'] in (g_dict[i] + "_sub_N", g_dict[i] + "_sub_P"):
                new_seq += gene_tables[i].iloc[j]['sequence']
                new_seq += new_second_part
            #print len(new_seq)
            concat[current_gene].append(new_seq)
        # print len(concat[current_gene])
    for gen in concat:

        #writer = csv.writer(f, delimiter= '\t', lineterminator='\n')
        for item in concat[gen]:
            f.write(gen + '\t' + item + '\n')
    f.close()
    pickle.dump(concat, open("concat.csv", "wb"))
    

sampling(mutated_genes)
