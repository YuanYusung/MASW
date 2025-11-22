'''
Multichannel Analysis of Surface Waves (MASW)
Modified based on Hongrui Qiu's codes by Yusong Yuan on Sep 4 18:00 2025 
'''
import numpy as np
from scipy.fftpack import fft,fftfreq# libraries for FFT analysis
import os
import glob
from obspy import read
import pickle
import matplotlib.pyplot as plt

def MASW(x,t,fs,datain,ct=np.linspace(0.2,2.5,231),
         fmin=0.5,fmax=5):
    '''
    Multi-channel Analysis of Surface Wave
    
    Input:
        x - distance, 1D array
        t - time axis
        fs - sampling rate, Hz
        datain - input wavefield , 2D array
            waveform at a certain receiver: datain[i,:]
            wavefield snapshot at a certain lapse time: datain[:,j]
        ct - grid search phase velocity, 1D array, km/s
        fmin,fmax - min and max of the frequency range for MASW
        
    Output:
        f - output frequency axis, 1D array, for plotting
        AsSum - MASW diagram, 2D array
            intensity distribution at certain frequency: AsSum[i,:]
            intensity distribution at certain grid search velocity: AsSum[:,j]
            
    This code is modified based on the Matlab script here:
        https://www.mathworks.com/matlabcentral/fileexchange/65138-masw-dispersion-curve
    '''
    # zero padding
    npts = datain.shape[1]; zeropad = np.zeros([datain.shape[0],5*npts])
    zeropad[:,0:npts] = datain
    # get the frequency axis
    freq = fftfreq(5*npts,d=1./fs); indf = (freq > fmin) & (freq < fmax)
    f = freq[indf]
    
    Ufft = fft(zeropad,axis=-1)[:,indf]
    w = f*2*np.pi
    
    Rnorm = Ufft/abs(Ufft)
    AsSum = np.zeros([len(w),len(ct)])
    for ii in range(len(w)):
        if ii % 100 == 0:
            print (ii,len(w))
        As = np.zeros([len(x),len(ct)],dtype=complex)
        for n in range(len(x)):
            As[n,:] = np.exp(1j*w[ii]*(x[n]/ct))*Rnorm[n,ii]
        AsSum[ii,:] = abs(np.sum(As,axis=0))
        AsSum[ii,:] /= max(AsSum[ii,:])
        
    return ct,f,AsSum

def read_sac():
    fout = 'All_sacdata_ZZ.pickle'
    vmax = 5.; vmin = 0.2
    if os.path.isfile(fout):
        st = np.load(fout,allow_pickle=True)
        return st
    dirpath = 'Data_CCFs/'
    sac_files = glob.glob(os.path.join(dirpath, "*.sac"))  

    st = None
    for fin in sac_files:
        if st is None:
            st = read(fin)
        else:
            st += read(fin)

    for i, tr in enumerate(st):
        tr.data = (tr.data + tr.data[::-1]) / 2 
        dist = tr.stats.sac.dist
        tminf = dist/vmax; tmaxf = dist/vmin
        st[i] = tapering(tr, tminf, tmaxf)

    st.filter('bandpass',freqmin=0.3,freqmax=3,zerophase=True)
    with open(fout, 'wb') as handle:
        pickle.dump(st, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return st

def tapering(tr,tmin,tmax):
    '''
    subroutine for tapering utilizing obspy
    
    tr - obspy trace object
    taxis - time axis of the input trace
    tmin,tmax - the start and end of the tapering window
    '''
    npts = tr.stats.sac.npts; t0 = tr.stats.sac.b; delta = tr.stats.delta
    taxis = np.linspace(t0,t0+npts*delta-delta,npts)
    ind_taper = (taxis >= tmin) & (taxis <= tmax)
    
    tr_copy = tr.copy()
    tr_copy.data = tr.data[ind_taper]
    tr_copy.taper(max_percentage=0.05)
    
    data = np.zeros_like(tr.data)
    data[ind_taper] = tr_copy.data
    tr_copy.data = data
    return tr_copy

def MASW_avg_dispersion():

    st = read_sac()
    wfs=[]; dist_lst=[]
    for tr in st:
        dist = tr.stats.sac.dist
        wf = tr.data / np.max(np.abs(tr.data))
        wfs.append(wf)
        dist_lst.append(dist)

    tr = st[0]; t0 = tr.stats.sac.b
    npts = tr.stats.npts; delta = tr.stats.delta
    taxis = np.linspace(t0,t0+npts*delta-delta,npts)
    
    fs = tr.stats.sampling_rate

    wfs_arr = np.array(wfs); dist_arr = np.array(dist_lst)
    indsort = np.argsort(dist_arr)

    fname = 'MASW_raw_data'
    if os.path.isfile(fname+'.npz'):
                out = np.load(fname+'.npz')
                ct = out['ct']; f = out['f']; AsSum = out['AsSum']
    else:
        ind = dist_arr[indsort] > 0.6 # exclude station pairs that are too close
        ct,f,AsSum = MASW(dist_arr[indsort][ind],taxis,fs,wfs_arr[indsort[ind],:],
                    ct=np.linspace(2.0,4.5,231),fmin=0.3,fmax=4.)
        out = {'ct':ct,'f':f,'AsSum':AsSum}
        np.savez(fname,**out)
        
    return ct,f,AsSum

def plot_cc():
    st = read_sac()
    st.filter('bandpass',freqmin=0.3,freqmax=3,zerophase=True)

    wfs, dist_lst = [], []
    for tr in st:
        dist = tr.stats.sac.dist
        wf = tr.data / np.max(np.abs(tr.data))
        wfs.append(wf)
        dist_lst.append(dist)

    wfs = np.array(wfs)
    dist_arr = np.array(dist_lst)

    tr0 = st[0]
    t0 = tr0.stats.sac.b
    npts = tr0.stats.npts
    delta = tr0.stats.delta
    taxis = np.linspace(t0, t0 + (npts - 1) * delta, npts)

    ind = np.argsort(dist_arr)
    X, Y = np.meshgrid(taxis, dist_arr[ind])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pcolormesh(X, Y, wfs[ind, :], cmap='coolwarm',
                    rasterized=True, shading='auto', vmin=-0.8, vmax=0.8)

    ymin, ymax = ax.get_ylim()
    t_line = np.array([ymin/3.0, ymax/3.0])
    ax.plot(t_line, [ymin, ymax], 'k--', lw=1.2)
    ax.annotate('3 km/s', xy=(1.5, ymin*0.8 + ymax*0.2),
                rotation=63, color='k', fontsize=18)

    ax.set_xlim(-1, 8)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('Correlation time (s)', fontsize=18)
    ax.set_ylabel('Interstation distance (km)', fontsize=18)
    ax.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig('Figure1_wavefield.png', dpi=300)

def plot_avg_disp():
    ct,f,AsSumR = MASW_avg_dispersion()
    X,Y = np.meshgrid(ct,f)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pcolormesh(Y,X,AsSumR,cmap='coolwarm',shading='auto',rasterized=True)
    ax.set_xlim(0.3,3)
    ax.set_xlabel('Freqency (Hz)',fontsize=18)
    ax.set_ylabel('Phase Velocity (km/s)',fontsize=18)
    fig.tight_layout()
    fig.savefig('Figure2_masw_avg_dispersion.png', dpi=300)

plot_cc()
plot_avg_disp()