# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 23:58:53 2017
"""

import numpy as np
import re

class Gene_transform:
    """Transform genes into Voss, z-curve, tetrahedron.
        Input for functions is string of {a, t, g, c} in lower case."""
    
    def __init__(self):
        self.NOT_GENE_CODE = re.compile(r'[^acgt]')
    
    def check_is_gene(self, genes):
        assert type(genes) is str
        return not bool(self.NOT_GENE_CODE.search(genes))
    
    """for string of lenght N returns array of shape (4, N)"""
    
    def voss(self, genes):
        assert type(genes) is str
        x = np.array(list(genes))
        u_a = (x == 'A').astype(int)
        u_t = (x == 'C').astype(int)
        u_g = (x == 'G').astype(int)
        u_c = (x == 'T').astype(int)
        return np.array([u_a, u_t, u_g, u_c])

    """is voss - flag for precomputed voss; if true => genes is voss array else genes is string;
       for string of lenght N returns array of shape (3, N)"""
    
    def z_curve(self, genes, is_voss=False, i = 1):
        if(not is_voss):
            voss = self.voss(genes)
        else:
            voss = genes
        assert type(voss) is np.ndarray
        assert voss.shape[0] == 4

        x_r = voss[0] - voss[1] + voss[2] - voss[3]
        x_g = voss[0] + voss[1] - voss[2] - voss[3]
        x_b = voss[0] - voss[1] - voss[2] + voss[3]
        
        res = np.array([x_r, x_g, x_b])  
        res.tofile("obs/" + str(i) + ".txt",",","%s") #export for txt
        
        return res
    
    def tetrahedron(self, genes, is_voss=False):
        if(not is_voss):
            voss = self.voss(genes)
        else:
            voss = genes
        assert type(voss) is np.ndarray
        assert voss.shape[0] == 4

        x_r = 2.0**0.5/3.0 * (2.0 * voss[3] - voss[1] - voss[2])
        x_g = 6.0**0.5/3.0 * (voss[1] - voss[2])
        x_b = 1.0/3.0 * (3.0 * voss[0] - voss[1] - voss[2] - voss[3])
        return np.array([x_r, x_g, x_b])  
    
    def z_curve_to_img(self, z):
        assert type(z) is np.ndarray
        assert z.shape[0] == 3
        image = z.T
        height = int(image.shape[0]**0.5)
        image = image[:height**2].reshape(height, height, 3)
        image = ((image + 1) * 255/2).astype("uint8")
        return image

    def tetrahedron_to_img(self, z):
        return self.z_curve_to_img(z)

        
        