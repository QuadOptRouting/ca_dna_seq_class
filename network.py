# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 00:16:34 2017

@author: 
"""
import torch
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
from torch.autograd import Variable
from torch import nn
import numpy as np
import matplotlib.pyplot as plt
import os.path
from skimage import io, transform

class GenomesDataset(Dataset):
    def __init__(self, img_arr):
        '''img_arr is an array of images
            assuming they all have the same shape
            and structure channel-width-heigth'''
        self.img_arr = img_arr
        self.ch, self.w, self.h = img_arr[0].shape
    
    def __len__(self):
        return len(self.img_arr)
    
    def __getitem__(self, idx):
        return self.img_arr[idx]

class ProjectNetwork(nn.Module):
    def __init__(self, input_size, input_channels):
        super(ProjectNetwork, self).__init__()
        self.conv_layers = nn.Sequential(nn.Conv2d(input_channels, 8, (3, 3), padding=1),
                                         nn.ReLU(),
                                         nn.MaxPool2d(2),
                                         nn.Conv2d(8, 16, (3, 3), padding=1),
                                         nn.ReLU(),
                                         nn.MaxPool2d(2)
                                        )
        self.linear_layers = nn.Sequential(nn.Linear(input_size * input_size / 4 * 16, 10),
                                           nn.LogSoftmax())
    
    def forward(self, x):
        x = self.conv_layers(x)
        x = x.view(x.size(0), -1)
        x = self.linear_layers(x)
        return x

