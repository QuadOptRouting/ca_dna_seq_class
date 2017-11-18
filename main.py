# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 00:10:23 2017

@author:
"""
import network
def main():
    X = []
    i = 1
    while os.path.isfile("obs/" + str(i) + ".txt"):
        if i >= 1000000:
            break

        observ = np.loadtxt("obs/" + str(i) + ".txt", delimiter = ",", dtype = np.int)
        observ = observ.reshape(3,len(observ)/3)

        X.append(observ)
        i += 1
    
    dataset = GenomesDataset(X)
    model = ProjectNetwork(...).cuda()
    loss_fn = torch.nn.NLLLoss().cuda()
    
if __name__ == '__main__':
    main()
