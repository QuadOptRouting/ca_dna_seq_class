import numpy as np
import re
import random
import matplotlib.pyplot as plt


class Gene_transform:
    """Transform genes into Voss, z-curve, tetrahedron.
        Input for functions is string of {a, t, g, c} in lower case."""
    
    def __init__(self):
        self.NOT_GENE_CODE = re.compile(r'[^acgt]')
    
    def check_is_gene(self, genes):
        assert isinstance(genes, str)
        return not bool(self.NOT_GENE_CODE.search(genes))
    
    """for string of lenght N returns array of shape (4, N)"""
    
    def voss(self, genes):
        assert self.check_is_gene(genes)
        x = np.array(list(genes))
        u_a = (x == 'a').astype(int)
        u_t = (x == 'c').astype(int)
        u_g = (x == 'g').astype(int)
        u_c = (x == 't').astype(int)
        return np.array([u_a, u_t, u_g, u_c])

    """is voss - flag for precomputed voss; if true => genes is voss array else genes is string;
       for string of lenght N returns array of shape (3, N)"""
    
    def z_curve(self, genes, is_voss=False):
        if(not is_voss):
            voss = self.voss(genes)
        else:
            voss = genes
        assert isinstance(voss, np.ndarray)
        assert voss.shape[0] == 4
        x_r = voss[0] - voss[1] + voss[2] - voss[3]
        x_g = voss[0] + voss[1] - voss[2] - voss[3]
        x_b = voss[0] - voss[1] - voss[2] + voss[3]
        return np.array([x_r, x_b, x_g])
    
    def tetrahedron(self, genes, is_voss=False):
        if(not is_voss):
            voss = self.voss(genes)
        else:
            voss = genes
        assert isinstance(voss, np.ndarray)
        assert voss.shape[0] == 4
        x_r = 2.0**0.5/3.0 * (2.0 * voss[3] - voss[1] - voss[2])
        x_g = 6.0**0.5/3.0 * (voss[1] - voss[2])
        x_b = 1.0/3.0 * (3.0 * voss[0] - voss[1] - voss[2] - voss[3])
        return np.array([x_r, x_b, x_g])  
    
    """cases: small_square, big_square, rectangle"""
    def z_curve_to_img(self, z, case="small_square"):
        assert isinstance(z, np.ndarray)
        assert z.shape[0] == 3
        image = z.T
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
        image = ((image + 1) * 255/2).astype("uint8")
        return image

    def tetrahedron_to_img(self, z):
        return self.z_curve_to_img(z)
        
        
class Gene_transform_arr:
    """Transform array of genes into Voss, z-curve, tetrahedron.
        Input for functions is array of strings of {a, t, g, c} in lower case."""
    
    def __init__(self):
        self.NOT_GENE_CODE = re.compile(r'[^acgt]')
    
    def check_is_gene(self, genes):
        assert isinstance(genes, np.ndarray)
        for string in genes:
            if(bool(self.NOT_GENE_CODE.search(string))):
                return False
        return True
    
    """for array of strings with shape (M, N) returns array with shape (4, M, N)"""
    
    def voss(self, genes):
        assert self.check_is_gene(genes)
        f = lambda x: list(x)
        x = np.array(map(f, genes))    
        u_a = (x == 'a').astype(int)
        u_t = (x == 'c').astype(int)
        u_g = (x == 'g').astype(int)
        u_c = (x == 't').astype(int)
        return np.array([u_a, u_t, u_g, u_c])

    """is voss - flag for precomputed voss; if true => genes is voss array else genes is string;
       for array of strings with shape (M, N) returns array with shape (3, M, N)"""
    
    def z_curve(self, genes, is_voss=False):
        if(not is_voss):
            voss = self.voss(genes)
        else:
            voss = genes
        assert isinstance(genes, np.ndarray)
        assert voss.shape[0] == 4
        x_r = voss[0] - voss[1] + voss[2] - voss[3]
        x_g = voss[0] + voss[1] - voss[2] - voss[3]
        x_b = voss[0] - voss[1] - voss[2] + voss[3]
        return np.array([x_r, x_b, x_g])
    
    def tetrahedron(self, genes, is_voss=False):
        if(not is_voss):
            voss = self.voss(genes)
        else:
            voss = genes
        assert isinstance(voss, np.ndarray)
        assert voss.shape[0] == 4
        x_r = 2.0**0.5/3.0 * (2.0 * voss[3] - voss[1] - voss[2])
        x_g = 6.0**0.5/3.0 * (voss[1] - voss[2])
        x_b = 1.0/3.0 * (3.0 * voss[0] - voss[1] - voss[2] - voss[3])
        return np.array([x_r, x_b, x_g])
    
    """cases: small_square, big_square, rectangle"""
    def z_curve_to_img(self, z, case = "small_square"):
        assert isinstance(z, np.ndarray)
        assert z.shape[0] == 3
        Array_list = np.split(z, z.shape[1], axis=1).reshape((z.shape[0], z.shape[2]))
        Image_list = []
        for image in Array_list:
            image = image.T
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
            image = ((image + 1) * 255/2).astype("uint8")
            Image_list.append(image)
        return Image_list

    def tetrahedron_to_img(self, z):
        return z_curve_to_img(z)
        
        

def usage_example():
	lenght = 1000000
	genes = ""
	for _ in xrange(lenght):
		seed = random.randint(1, 4)
		if(seed == 1):
			genes += "a"
		elif(seed == 2):
			genes += "t"
		elif(seed == 3):
			genes += "g"
		else:
			genes += "c"

	print len(genes)
	print genes[200:240]

	transformer = Gene_transform()
	print transformer.check_is_gene(genes)
	voss = transformer.voss(genes)
	z_curve = transformer.z_curve(voss, is_voss=True)
	tetrahedron = transformer.tetrahedron(voss, is_voss=True)
	img = transformer.z_curve_to_img(z_curve)
	print img.shape
	plt.imshow(img)
	plt.show()


if __name__ == "__main__":
	usage_example()      

