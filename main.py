# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 00:10:23 2017

@author:
"""
from torch.utils.data import Dataset, DataLoader
from torch.autograd import Variable
import network
import torch
import os
import numpy as np
import tqdm


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
    
    dataset = network.GenomesDataset(X)
    data_loader = DataLoader(dataset, batch_size=32, shuffle=True)
    model = network.ProjectNetwork().cuda()
    loss_fn = torch.nn.NLLLoss().cuda()
    learning_rate = 1e-4
    optimizer = torch.optim.Adam(params=model.parameters(), lr=learning_rate)
    nepoch = 10
    tqt = tqdm(xrange(nepoch), desc="Train")

    all_loss = {}
    test_loss = {}

    for t in tqt:
        print t
        j = 0
        ep_loss = 0
        for X, y in data_loader:
            X = Variable(X).type(torch.cuda.IntTensor)
            y = Variable(y).type(torch.cuda.IntTensor)  ## for classification - vector of longs

            y_pred = model(X)
            loss = loss_fn(y_pred, y)
            ep_loss += loss.data[0]

            optimizer.zero_grad()

            loss.backward()

            optimizer.step()
            j += 1
        all_loss.append(ep_loss / j)
        tqt.set_description("Loss = {}".format(all_loss[-1]))


if __name__ == '__main__':
    main()
