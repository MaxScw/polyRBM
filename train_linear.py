import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.autograd import Variable
import time
import os
from tqdm import tqdm

import polyRBM
import analysisUtility as anaU


# training
chainLengths = [8, 32]
networkSizes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]
parent = os.getcwd()
savepath = ###directory to save training results in###
fname_8 = ###file path of N=8 simulation data###
fname_32 = ###file path of N=32 simulation data###
fnames = [fname_8, fname_32]

for ind, l in enumerate(chainLengths):
    fname = fnames[ind]
    dataset = polyRBM.CubicChainDataset(fname, 2)
    div = int(0.7*len(dataset.data))
    testset = dataset.bits[div:]
    dataset.data = torch.from_numpy(dataset.bits[:div]).float()
    
    for n in tqdm(networkSizes):
        run = f'rw_exclVol_N{l}_h{int(n)}' #unique identifier is automatically created
        os.mkdir(parent+'/'+savepath+f'repr_results[model={run}]')
        batchsize = 10
        n_vis = 2*(l-2)
        n_hidden = int(n)
        k = 2
        rbm_lr = 10**-2
        rbm_epochs= 200
        rbm_verbose = False
        
        train_loader = polyRBM.DataLoader(dataset, batch_size=batchsize, shuffle=True) 
        rbm = polyRBM.RBM(n_vis=n_vis, n_hin=n_hidden, k=k, oneHot=False, chainlength=l)
        
        print(f"START training of RBM [N={l}, n_hidden={int(n)}]")
        print(os.getcwd())
        loss, recon_error, countList = polyRBM.train_RBM(rbm=rbm, train_loader=train_loader, epochs=rbm_epochs, 
                                                 lr=rbm_lr, verbose=rbm_verbose, ninputs=n_vis, 
                                                 save_model=True, save_path=savepath+f'{run}_rbm', 
                                                 testset=torch.from_numpy(testset).float(), momentum=0.5, calc_loops=True)
        np.save(savepath+f'repr_results[model={run}]/{run}_loss.npy', loss)
        np.save(savepath+f'repr_results[model={run}]/{run}_loops.npy', countList)
        all_samples = []
        
        for param in rbm.parameters():
            param.requires_grad = False
        
        for i in tqdm(range(10)):
            initialize = torch.from_numpy(np.random.binomial(1, 0.5, (1000, (l-2)*2))).float()
            rbm.k = 800
            samples = rbm(initialize)[2].detach().numpy()
            all_samples.append(samples)
        samples = np.concatenate(all_samples)
        np.save(savepath+f'repr_results[model={run}]/{run}_samples.npy', samples)
