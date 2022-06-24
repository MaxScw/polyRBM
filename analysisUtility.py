import numpy as np
import polyRBM 
from matplotlib import pyplot as plt


# function to load files according to savepath + unique identifier

def loader(savepath, run, rbm=False, dataset=False, data=False, dataTF=False, 
           samples=False, samplesTF=False, samplesTFsorted=False, loss=False,
           dataTFVEC=False, loops=False):
    path = savepath + 'repr_results[model={}]/{}_'.format(run, run)
    if rbm==True:
        rbm = polyRBM.torch.load(savepath+'{}_rbm'.format(run))
        return rbm
    if dataset==True:
        dataset = polyRBM.torch.load(savepath+'{}_dat'.format(run))
        return dataset
    if data==True:
        dpath = path + 'data.npy'
        data = np.load(dpath, allow_pickle=True)
        return data
    if dataTF==True:
        dtpath = path + 'data_tf.npy'
        dataTF = np.load(dtpath, allow_pickle=True)
        return dataTF
    if dataTFVEC==True:
        dtvpath = path + 'data_tf_vec.npy'
        dataTFVEC = np.load(dtvpath, allow_pickle=True)
        return dataTFVEC
    if samples==True:
        spath = path + 'samples.npy'
        samples = np.load(spath)
        return samples
    if samplesTF==True:
        stpath = path + 'samples_tf.npy'
        samplesTF = np.load(stpath, allow_pickle=True)
        return samplesTF
    if samplesTFsorted==True:
        stspath = path + 'samples_tf_sorted.npy'
        samplesTFsorted = np.load(stspath, allow_pickle=True)
        return samplesTFsorted
    if loss==True:
        lpath = path + 'loss.npy'
        print(lpath)
        loss = np.load(lpath)
        return loss
    if loops==True:
        loopath = path + 'loops.npy'
        print(loopath)
        loops = np.load(loopath)
        return loops


# function to calculate number of loops in dataset

def loop_count(data):
    size = len(data)
    samples = polyRBM.bitsToBonds2DYu(data)
    end_pos = np.cumsum(samples, axis=1)[:, -1, :]
    tmploopcounts = np.zeros((1000, 100))            
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
                
    loopcount = np.sum(tmploopcounts[:, :], axis=0)/size
    return loopcount


# function to calculate the squared radius of gyration
# on a dataset for nsamples samples

def rgx_calc(data, nsamples):
    copy = data[:nsamples]
    copy = np.cumsum(copy, axis=1)
    rows = len(copy[:, 0])
    cols = len(copy[0, :])
    
    rgx = []
    rgy = []
    rgz = []
    
    for i in range(rows):
        mean = np.mean(1.*copy[i, :])
        
        tmp_rgx = []
        tmp_rgy = []
        tmp_rgz = []
        for j in range(cols):
            tmp_rgx.append((copy[i, j][0] - mean[0])**2)
            tmp_rgy.append((copy[i, j][1] - mean[1])**2)
            tmp_rgz.append((copy[i, j][2] - mean[2])**2)
        rgx.append(np.mean(np.array(tmp_rgx)))
        rgy.append(np.mean(np.array(tmp_rgy)))
        rgz.append(np.mean(np.array(tmp_rgz)))
       
    rgx = np.array(rgx)
    rgy = np.array(rgy)
    rgz = np.array(rgz)
    rgtot = np.sum((rgx, rgy, rgz), axis=0)
    return  rgx, rgy, rgz, rgtot


# function to calculate the squared end-to-end distance
# for a given dataset and nsamples samples

def endToend_calc(data, nsamples):
    copy = data[:nsamples]
    end = np.sum(data, axis=1)
    endToend = np.zeros(len(end))
    for i in range(len(end)):
        endToend[i] = np.sum(np.sum(end[i], axis=0)**2)
    return endToend


# function to calculate bond vector correlation of dataset for
# nsamples samples
    
def bondvec_correlation(data, nsamples):
    data = data[:nsamples]
    rows = len(data[:, 0])
    cols = len(data[0, :])
    b_corr = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            b_corr[i, j] = np.dot(data[i, 0], data[i, j])
    return np.mean(b_corr, axis=0)


# function for calculating the bond angle distribution for dataset
# and nsamples samples
  
def angles(data, nsamples):
    data = data[:nsamples]
    rows = len(data[:, 0])
    cols = len(data[0, :])-1
    scalar = np.zeros((rows, cols))
    
    count = 0
    
    multi = data[:, 1:]*data[:, :-1]
    mag1 = data[:, 1:]*data[:, 1:]
    mag2 = data[:, :-1]*data[:, :-1]
    for i in range(rows):
        for j in range(cols):
            scalar[i, j] = np.sum(multi[i, j])/(np.sqrt(np.sum(mag1[i, j]))*np.sqrt(np.sum(mag2[i, j])))
            if abs(scalar[i, j]-1) < 10**-10:
                scalar[i, j] = 0.999
            elif abs(scalar[i, j]+1) < 10**-10:
                count += 1
                scalar[i, j] = -0.999
            scalar[i, j] = np.arccos(scalar[i, j])
    
    return scalar, count/(rows*cols) 


# function for calculating number of conformations present in
# reconstruction data set (dataset) but not in original data set (dataset)

def not_in_set(dataset, recon):
    not_in_set = 0
    unD = np.unique(dataset, axis=0)
    unR, unR_count = np.unique(recon.astype(int), return_counts=True, axis=0)
    
    print(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) in unD )
    for ind, u in enumerate(unR):
        boo = (u==unD).all(axis=1)
        ind = np.where(boo==True)[0]
        if ind.size > 0:
            not_in_set += unR_count[ind[0]]
    return not_in_set/len(recon)
    

'''more complex utility, e.g. plot functions'''


# plot training loss, e.g. <v-h>_0 - <v-h>_k

def plot_loss(loss, start=0):
    fig, ax = plt.subplots(1)
    fig.set_size_inches(7, 4)

    ax.plot(loss[start:], label="<v$\cdot$p$_h$>$_{data}$ - <v$\cdot$p$_h$>$_{recon}$")
    ax.set_xlabel('epochs', fontsize='x-large')
    ax.legend(fontsize='x-large')
    

# plot comparison of radius of gyration distribution between
# data (data) and reconstruction (recon)
    
def plot_Rg2(data, recon):
    fig, ax = plt.subplots(4)
    fig.tight_layout()
    fig.set_size_inches(10, 16)
    means_dat = np.mean(data, axis=1)
    means_recon = np.mean(recon, axis=1)

    for i, iax in enumerate(ax):
        iax.set_xlabel("Rg$^2$")
        iax.set_ylabel("%")
        iax.hist(data[i], bins=100, range=(0, 4*means_dat[-1]), label='data', density=True)
        iax.hist(recon[i], bins=100, range=(0, 4*means_dat[-1]),  label='reconstructed', density=True, histtype='step')
        iax.axvline(x=means_dat[i], ls='dotted', c='blue', label='<Rg$^2$>$_{dat}$')
        iax.axvline(x=means_recon[i], ls='dotted', c='red', label='<Rg$^2$>$_{recon}$')
        iax.legend()
    return means_dat, means_recon
    
 
# plot comparison of distribution of squared end-to-end distance between
# data (data) and reconstruction (recon) 
    
def plot_Re2(data, recon):
    fig, ax = plt.subplots(1)
    fig.set_size_inches(10, 6)
    mean_dat = np.mean(data)
    mean_recon = np.mean(recon)

    ax.set_xlabel("R$_e^2$")
    ax.set_ylabel("%")
    ax.hist(data, bins=20, range=(0, 4*mean_recon), label='data', density=True)
    ax.hist(recon, bins=20, range=(0, 4*mean_recon), label='reconstructed', histtype='step', color='orange', density=True)
    ax.axvline(x=mean_dat, c='blue', ls='dotted', label='<R$_e^2$>$_{data}$')
    ax.axvline(x=mean_recon, c='red', ls='dotted', label='<R$_e^2$>$_{recon}$')
    ax.legend()
    
    
# plot comparison between bond-vector correlation of data (data)
# and reconstruction (samples) for nsamples conformations
    
def plot_corr(data, samples, nsamples=0):
    if nsamples==0:
        corr_data = bondvec_correlation(data, len(data[:, 0]))
        corr_samples = bondvec_correlation(samples, len(samples[:, 0]))
    else:
        corr_data = bondvec_correlation(data, nsamples)
        corr_samples = bondvec_correlation(samples, nsamples)
    fig, ax = plt.subplots(1)
    
    ax.plot(np.linspace(1, len(corr_data), len(corr_data)), corr_data, label='data')
    ax.plot(np.linspace(1, len(corr_samples), len(corr_samples)), corr_samples, label='reproduction')
    ax.scatter(np.linspace(1, len(corr_data), len(corr_data)), corr_data, label='data')
    ax.scatter(np.linspace(1, len(corr_samples), len(corr_samples)),corr_samples, label='reproduction')
    ax.set_xlabel("i")
    ax.set_ylabel("<$\\vec{b_0}\\cdot\\vec{b_i}$>")
    ax.set_xlim(0.5, len(corr_data)+0.5)
    fig.legend()
    fig.set_size_inches(10, 8)
    fig.tight_layout()

    
# plot comparison between bond angle distribution between data (data)
# and reconstruction (recon)
    
def plot_angleDistribution(data, recon, nsamples=0, scale=1, ymax=0):
    if nsamples==0:
        angles_d, count_d = angles(data, len(data))
        angles_r, count_r = angles(recon, len(recon))
    else:
        angles_d, count_d = angles(data, nsamples)
        angles_r, count_r = angles(recon, nsamples)
    fig, ax = plt.subplots(1)
    ax.hist(angles_d.flatten()*180/np.pi, bins=200, range=(-5, 185), label='data', density=True)
    ax.hist(angles_r.flatten()*180/np.pi, bins=200, range=(-5, 185), histtype='step', label='reconstructed', density=True)
    ax.set_xticks(np.arange(0, 180+1, 10))
    fig.legend(fontsize=12*scale)
    fig.set_size_inches(10*scale, 5*scale)
    if ymax!=0:
        ax.set_ylim(0, ymax)
    return count_d, count_r
