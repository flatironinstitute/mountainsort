import numpy as np
from synthesize_single_waveform import synthesize_single_waveform

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import readmda,writemda32

processor_name='pyms.synthesize_random_waveforms'
processor_version='0.1'
def synthesize_random_waveforms(*,waveforms_out=None,geometry_out=None,M=5,T=500,K=20,upsamplefac=13):
    """
    Synthesize random waveforms for use in creating a synthetic timeseries dataset

    Parameters
    ----------
    waveforms_out : OUTPUT
        Path to waveforms mda file. Mx(T*upsamplefac)xK
    geometry_out : OUTPUT
        (Optional) Path to geometry csv file
    M : int
        (Optional) Number of channels
    T : int
        (Optional) Number of timepoints for a waveform, before upsampling
    K : int
        (Optional) Number of waveforms to synthesize
    upsamplefac : int
        (Optional) used for upsampling the waveforms to avoid discretization artifacts
    """    
    geometry=None
    avg_durations=[200,10,30,200]
    avg_amps=[0.5,10,-1,0]
    rand_durations_stdev=[10,4,6,20]
    rand_amps_stdev=[0.2,3,0.5,0]
    rand_amp_factor_range=[0.5,1]
    geom_spread_coef1=0.2
    geom_spread_coef2=1
    average_peak_amplitude=10
    timeshift_factor=3
    
    if not geometry:
        geometry=np.zeros((2,M))
        geometry[0,:]=np.arange(1,M+1)
    
    geometry=np.array(geometry)
    avg_durations=np.array(avg_durations)
    avg_amps=np.array(avg_amps)
    rand_durations_stdev=np.array(rand_durations_stdev)
    rand_amps_stdev=np.array(rand_amps_stdev)
    rand_amp_factor_range=np.array(rand_amp_factor_range)
    
    neuron_locations=get_default_neuron_locations(M,K,geometry)

    ## The waveforms_out
    WW=np.zeros((M,T*upsamplefac,K))
    
    for k in range(1,K+1):
        for m in range(1,M+1):
            diff=neuron_locations[:,k-1]-geometry[:,m-1]
            dist=np.sqrt(np.sum(diff**2))
            durations0=np.maximum(np.ones(avg_durations.shape),avg_durations+np.random.randn(1,4)*rand_durations_stdev)*upsamplefac
            amps0=avg_amps+np.random.randn(1,4)*rand_amps_stdev
            waveform0=synthesize_single_waveform(N=T*upsamplefac,durations=durations0,amps=amps0)
            waveform0=np.roll(waveform0,int(timeshift_factor*dist*upsamplefac))
            waveform0=waveform0*np.random.uniform(rand_amp_factor_range[0],rand_amp_factor_range[1])
            WW[m-1,:,k-1]=waveform0/(geom_spread_coef1+dist*geom_spread_coef2)

    peaks=np.max(np.abs(WW),axis=(0,1))
    WW=WW/np.mean(peaks)*average_peak_amplitude

    if waveforms_out:
        writemda32(WW,waveforms_out)
        if geometry_out:
            np.savetxt(geometry_out,geometry.transpose(),delimiter=",",fmt="%g")
            return True
        else:
            return True
    else:
        return (WW,geometry)
    
def get_default_neuron_locations(M,K,geometry):
    num_dims=geometry.shape[0]
    neuron_locations=np.zeros((num_dims,K))
    for k in range(1,K+1):
        if K>0:
            ind=(k-1)/(K-1)*(M-1)+1
            ind0=int(ind)
            if ind0==M:
                ind0=M-1
                p=1
            else:
                p=ind-ind0
            if M>0:
                neuron_locations[:,k-1]=(1-p)*geometry[:,ind0-1]+p*geometry[:,ind0]
            else:
                neuron_locations[:,k-1]=geometry[:,0]
        else:
            neuron_locations[:,k-1]=geometry[:,0]
        
    return neuron_locations
    
def test_synthesize_random_waveforms():
    M,T,K = 5,800,20
    upsamplefac=13
    waveforms,geometry = synthesize_random_waveforms(M=M,T=T,K=K,upsamplefac=upsamplefac)
    print(waveforms.shape)
    #waveforms=readmda('tmp.waveforms.mda')
    #assert(waveforms.shape==(M,T*upsamplefac,K))
    #import matplotlib.pyplot as plt
    #for k in range(1,K+1):
    #    plt.figure()
    #    for m in range(1,M+1):
    #        plt.plot(waveforms[m-1,:,k-1])
    
    return True

synthesize_random_waveforms.test=test_synthesize_random_waveforms
synthesize_random_waveforms.name = processor_name
synthesize_random_waveforms.version = processor_version

if __name__ == '__main__':
    print ('Running test')
    test_synthesize_random_waveforms()
    
