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
from tqdm import tqdm

  
# PyTorch module containing the RBM infrastructure 
       
class RBM(nn.Module):
    def __init__(self, n_vis=10, n_hin=5, k=1, chainlength=0, momentum=0, regular=0):
        super(RBM, self).__init__()
        
        initial = np.sqrt(6/(n_vis + n_hin))
        self.W = nn.Parameter(torch.randn(n_hin, n_vis)*initial -initial/2)
        self.v_bias = nn.Parameter(torch.zeros(n_vis))
        self.h_bias = nn.Parameter(torch.zeros(n_hin))
        self.k = k
        self.sample_mode = False
        
        self.momentum = momentum
        self.regular = regular
        
        self.dW = 0
        self.da = 0
        self.db = 0
        
    def sample_from_p(self, p):
        return F.relu(torch.sign(p - Variable(torch.rand(p.size()))))
    
    def v_to_h(self, v):
        p_h = torch.sigmoid(F.linear(v, self.W, self.h_bias))
        sample_h = self.sample_from_p(p_h)
        return p_h,sample_h
    
    def h_to_v(self, h):
        p_v = torch.sigmoid(F.linear(h, self.W.t(), self.v_bias))
        sample_v = self.sample_from_p(p_v)
        return p_v, sample_v
    
    def prob_to_ind(self, p_sum):
        return torch.argmax(p_sum[p_sum<(torch.rand(1)*(1-0.0001)+0.0001)]).detach().numpy()
        
    def update(self, dW, da, db, lr):
        with torch.no_grad():
            self.dW *= self.momentum
            self.da *= self.momentum
            self.db *= self.momentum
            
            self.dW += dW*lr
            self.da += da*lr
            self.db += db*lr
            
            self.W += self.dW
            self.v_bias += self.da
            self.h_bias += self.db
            
            self.W.grad = None
            self.v_bias.grad = None
            self.h_bias.grad = None
            
            self.W = nn.Parameter(torch.clamp(self.W, -10, 10))
            
    def forward(self, v):
        pre_h1, h1 = self.v_to_h(v)
        h_ = h1
        for _ in range(self.k):
            pre_v_,v_ = self.h_to_v(h_)
            pre_h_,h_ = self.v_to_h(v_)
        return  v, h1, v_, h_, pre_h1, pre_h_
    
    def free_energy(self, v):
        vbias_term = v.mv(self.v_bias)
        wx_b = F.linear(v, self.W, self.h_bias)
        hidden_term = wx_b.exp().add(1).log().sum(1)
        return (-hidden_term - vbias_term).mean()


# training function implementing the approximate Likelihood ascent

def train_RBM(rbm, train_loader, epochs, lr, verbose, ninputs, save_model, save_path='', testset=None, calc_loops=False, momentum=0):
    global_loss = []
    global_meansq = []
    loopcounts = np.zeros((epochs, 100))
    eToe = np.zeros(epochs)
    tuned_down1 = False
    tuned_down2 = False
    rbm.momentum = momentum
    
    for epoch in tqdm(range(epochs)):
        t1 = time.time()
        tmp_v = []
        tmp_vk = []
        tmp_ph = []
        tmp_phk = []
        tmploopcounts = np.zeros((1000, 100))
        dWs = []
        
        
        for ind, data in enumerate(train_loader):
            v, h, v_k, h_k, p_h, p_h_k = rbm(data)
            tmp_v.append(v)
            tmp_vk.append(v_k)
            tmp_ph.append(p_h)
            tmp_phk.append(p_h_k)
            
            dW = ((torch.mm(v.t(), p_h) - torch.mm(v_k.t(), p_h_k)).T)/len(v)
            dWs.append(dW.detach().numpy())
            da = torch.mean(v - v_k, 0)
            db = torch.mean(p_h - p_h_k, 0)
        
            if (epoch >= epochs/2) and (tuned_down1==False):
                    lr = lr*0.01
                    tuned_down1 = True
            if (epoch >= 3*epochs/4) and (tuned_down2==False):
                    lr = lr*0.01
                    tuned_down2 = True
        
            rbm.update(dW, da, db, lr)
        
        if (epoch+1)%1==0:
            if calc_loops == True:
                bits = rbm(testset)[2].detach().numpy()
                size = len(bits)
                samples = bitsToBonds2DYu(bits)
                end_pos = np.cumsum(samples, axis=1)[:, -1, :]
                endToend_distance = np.sqrt(np.sum(end_pos**2, axis=1))
                
                pos = np.zeros((len(samples[:, 0, 0]), len(samples[0, :, 0])+1, 3))
                pos[:, 1:, :] = np.cumsum(samples, axis=1)
                samples = pos
                already_zero = (samples==0).all(axis=2)
            
                for i in range(len(samples[0, :, 0])):
                    diff = np.array([samples[:, j, :] - samples[:, i, :] for j in range(len(samples[0, :, 0]))])
                    diff[i] = samples[:, i, :]
                    diff_zero = (diff==0).all(axis=2).T
                    new_zero = np.logical_and(diff_zero, np.logical_not(already_zero))
                    if i==0:
                        for j in range(i+1, len(samples[0, :, 0])):
                            zero = already_zero[:, j]
                            dist = j-i
                            tmploopcounts[i, dist] = len(zero[zero==True])
                    else:
                        for j in range(i+1, len(samples[0, :, 0])):
                            zero = diff_zero[:, j][new_zero[:, j]]
                            dist = j-i
                            tmploopcounts[i, dist] = len(zero==True)
                
                loopcounts[epoch, :] = np.sum(tmploopcounts[:, :], axis=0)/size
                eToe[epoch] = np.mean(endToend_distance)
            
        global_loss.append(np.mean(np.array(dWs)))
        
        t2 = time.time()-t1
        if verbose==True:
            print("loss for {} epoch: {} \n step in {} seconds".format(epoch+1, global_loss[-1], t2))
            print("end-to-end dist: {}".format(eToe[epoch]))
            
    if save_model == True:
        torch.save(rbm, save_path)
    return global_loss, eToe, loopcounts
    

'''transformation utility for random walk chains'''


# function for creating position vectors of a conformation's
# monomers from the relative bond vectors

def bondsToConf(b):
    
    N = b.shape[1]
    
    conf=np.zeros((b.shape[0],N+1,3)).astype(int)
    
    for i in range( N+1 ):
        conf[:,i,:]=np.sum(b[:,:i,:],axis=1)

    return conf
    

# functio to transform 2d-bond vectors (XxYx3 with the last of three
# positions empty) to Yu-inspired bit encoding of conformations (XxYx1)

def bondsToBits2DYu(b):
    
    bits=np.zeros((b.shape[0],b.shape[1]-1,2)).astype(int)
    
    scalar=b[:,1:,0]*b[:,:-1,0]+b[:,1:,1]*b[:,:-1,1]
    cross =b[:,1:,0]*b[:,:-1,1]-b[:,1:,1]*b[:,:-1,0]
    
    bits[:,:,0]=(scalar==-1) + (cross == -1)
    bits[:,:,1]=(scalar==-1) + (cross ==  1)
    
    bits = bits.reshape((bits.shape[0],bits.shape[1]*bits.shape[2]))

    return bits.astype(int)


# inverse of 'bondsToBits2DYu'

def bitsToBonds2DYu(bits_):
    
    N     = bits_.shape[1]//2
    bits  = bits_.reshape((bits_.shape[0],N,2))
    
    bonds = np.zeros((bits_.shape[0],N+1,3)).astype(int)
        
    shape = (bits_.shape[0],N,1)
    
    back   = (1*(bits[:,:,0]==1) * (bits[:,:,1]==1)).reshape(shape)
    scalar = (1*(bits[:,:,0]==0) * (bits[:,:,1]==0)).reshape(shape)
    cross  = (1*((bits[:,:,0]==1) * (bits[:,:,1]==0))-1*((bits[:,:,0]==0) * (bits[:,:,1]==1))).reshape(shape)
    
    b0     = np.zeros((1,1,3))
    b0[0,0,1] = 1
    
    for i in np.arange(N+1):
        if ( i == 0 ):
            bonds[:,i,:]=b0
        else:
            
            b1   =bonds[:,i-1,:]
            bond =-1*back[:,i-1]*b1+scalar[:,i-1]*b1
   
            b1flip=b1[:,(1,0,2)]
            b1flip[:,0]*=-1
  
            bond += cross[:,i-1]*b1flip
            bonds[:,i,:]=bond
     
    return bonds


# function to transform XxYx3 matrix of conformations into
# XxY array of numpy vectors of dimension 3

def makeVecs(data):
    rows = len(data[:, 0, :])
    cols = len(data[0, :, :])
    dataTF = np.ndarray((rows, cols), dtype=object)
    for i in range(rows):
        for j in range(cols):
            dataTF[i, j] = data[i, j, :]
    return dataTF


# pytorch Dataset class to save training data in for evaluation
# by the RBM

class CubicChainDataset(Dataset):
    def __init__(self, filename, dim):

        t1    = time.time()

        self.dim = dim
        tmp   = pd.read_csv(filename,sep=",")
        tmp   = np.array(tmp.iloc[:,1:].values)

        N     = tmp.shape[1]//3
        self.bonds = tmp.reshape((tmp.shape[0],N,3))
        
        self.bits  = bondsToBits2DYu(self.bonds)

        self.data = torch.from_numpy(self.bits)
        self.data = self.data.type(torch.float32)

        t2 = time.time() - t1
        print("initialization complete after {} seconds".format(t2))

    def save_data(self):
        return 0

    def __len__(self):
        return len(self.data)

    def __getitem__(self,idx):
        return self.data[idx, :]
