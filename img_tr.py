import numpy as np
import re
import random
import matplotlib.pyplot as plt


class Gene_transform:
    """Transform genes into Voss, z-curve, tetrahedron.
        Input for functions is string of {A, T, G, C} in upper case."""
    
    def __init__(self, cumsum=True):
        self.NOT_GENE_CODE = re.compile(r'[^ACGT]')
        self.cumsum = cumsum
    
    """Use this function. gene_seq => image"""
    def transform_to_img(self, gen_seq, case="big_square"):
        z_curve = self._z_curve(gen_seq)
        return self._z_curve_to_img(z_curve, case=case)
    
    def _check_is_gene(self, genes):
        assert isinstance(genes, str)
        return not bool(self.NOT_GENE_CODE.search(genes))
    
    """for string of lenght N returns array of shape (4, N)"""
    
    def _voss(self, genes):
        assert self._check_is_gene(genes)
        x = ['x']
        x += list(genes)
        x = np.array(x)
        u_a = (x == 'A').astype(int)
        u_t = (x == 'C').astype(int)
        u_g = (x == 'G').astype(int)
        u_c = (x == 'T').astype(int)
        return np.array([u_a, u_t, u_g, u_c])

    """is voss - flag for precomputed voss; if true => genes is voss array else genes is string;
       for string of lenght N returns array of shape (3, N)"""
    
    def _z_curve(self, genes, is_voss=False, i = 1):
        if(not is_voss):
            voss = self._voss(genes)
        else:
            voss = genes
        assert isinstance(voss, np.ndarray)
        assert voss.shape[0] == 4
        if self.cumsum:
            voss = np.cumsum(voss, axis=1)
        x_r = voss[0] - voss[1] + voss[2] - voss[3]
        x_g = voss[0] + voss[1] - voss[2] - voss[3]
        x_b = voss[0] - voss[1] - voss[2] + voss[3]
        res = np.array([x_r, x_g, x_b])  
        res.tofile("obs/" + str(i) + ".txt",",","%s") #export for txt
        return res
    
    """Attention! No cumulative option in tetrahedron; It won't work properly, don't use it"""
    def _tetrahedron(self, genes, is_voss=False):
        if(not is_voss):
            voss = self._voss(genes)
        else:
            voss = genes
        assert isinstance(voss, np.ndarray)
        assert voss.shape[0] == 4
        x_r = 2.0**0.5/3.0 * (2.0 * voss[3] - voss[1] - voss[2])
        x_g = 6.0**0.5/3.0 * (voss[1] - voss[2])
        x_b = 1.0/3.0 * (3.0 * voss[0] - voss[1] - voss[2] - voss[3])
        return np.array([x_r, x_g, x_b])  
    
    """cases: small_square, big_square, rectangle"""
    def _z_curve_to_img(self, z_curve, case="big_square"):
        assert isinstance(z_curve, np.ndarray)
        assert z_curve.shape[0] == 3
        image = z_curve.T
        height = int(image.shape[0]**0.5)
        if(case == "big_square"):
            height += 1
            image = np.concatenate((image, -np.ones((2*height, 3))), axis=0)
            image = image[:height**2].reshape(height, height, 3)
        elif(case == "small_square"):
            image = image[:height**2].reshape(height, height, 3)
        elif(case == "rectangle"):
            print "not implemented yet"
            return None
        else:
            print "wrong args"
            return None
        if self.cumsum:
            seq_lenght = z_curve.shape[1]
        else:
            seq_lenght = 1
        image = ((image + seq_lenght) * 255/(2.0 * seq_lenght)).astype("uint8")
        return image
    def _tetrahedron_to_img(self, z):
        return self.z_curve_to_img(z)


def usage_example():
	lenght = 1000
	genes = ""
	for _ in xrange(lenght):
		seed = random.randint(1, 4)
		if(seed == 1):
			genes += "A"
		elif(seed == 2):
			genes += "T"
		elif(seed == 3):
			genes += "G"
		else:
			genes += "C"

	print len(genes)
	print genes[200:240]

	transformer = Gene_transform()
	img = transformer.transform_to_img(genes);
	print img.shape
	plt.imshow(img)
	plt.show()


if __name__ == "__main__":
	usage_example() 

