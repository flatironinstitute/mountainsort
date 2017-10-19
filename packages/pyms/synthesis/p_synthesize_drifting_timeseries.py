import numpy as np
import scipy.interpolate as spi

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from p_synthesize_timeseries import synthesize_timeseries

from mlpy import readmda,writemda32,writemda64

processor_name='pyms.synthesize_drifting_timeseries'
processor_version='0.1'
def synthesize_drifting_timeseries(*,
        firings,
        waveforms,
        timeseries_out=None,
        noise_level=1,
        samplerate=30000,
        duration=60,
        waveform_upsamplefac=1,
        amplitudes_row=0,
        num_interp_nodes=2
    ):
    """
    Synthesize a electrophysiology timeseries from a set of ground-truth firing events and waveforms, and simulating drift (linear for now)

    Parameters
    ----------
    firings : INPUT
        (Optional) The path of firing events file in .mda format. RxL where 
        R>=3 and L is the number of events. Second row is the timestamps, 
        third row is the integer labels
    waveforms : INPUT
        (Optional) The path of (possibly upsampled) waveforms file in .mda
        format. Mx(T*waveform_upsample_factor)*(2K), where M is the number of
        channels, T is the clip size, and K is the number of units. Each unit
        has a contiguous pair of waveforms (interpolates from first to second
        across the timeseries)
    
    timeseries_out : OUTPUT
        The output path for the new timeseries. MxN

    noise_level : double
        (Optional) Standard deviation of the simulated background noise added to the timeseries
    samplerate : double
        (Optional) Sample rate for the synthetic dataset in Hz
    duration : double
        (Optional) Duration of the synthetic dataset in seconds. The number of timepoints will be duration*samplerate
    waveform_upsamplefac : int
        (Optional) The upsampling factor corresponding to the input waveforms. (avoids digitization artifacts)
    amplitudes_row : int
        (Optional) If positive, this is the row in the firings arrays where the amplitude scale factors are found. Otherwise, use all 1's
    num_interp_nodes : int
        (Optional) For drift, the number of timepoints where we specify the waveform (Default 2)
    """
    
    if type(firings)==str:
        F=readmda(firings)
    else:
        F=firings

    if amplitudes_row==0:
        F=np.concatenate((F,np.ones((1,F.shape[1]))))
        amplitudes_row=F.shape[0]
        
    times=F[1,:]
    times_normalized=times/(duration*samplerate) #normalized between 0 and 1
    labels=F[2,:]
    amps=F[amplitudes_row-1,:]
            
    F=np.kron(F,[1]*num_interp_nodes) #duplicate every event!
    
    for j in range(num_interp_nodes):
        F[amplitudes_row-1,j::num_interp_nodes]=amps*time_basis_func(j,num_interp_nodes,times_normalized)
        # adjust the labels
        F[2,j::num_interp_nodes]=(labels-1)*num_interp_nodes+j+1 #remember that labels are 1-indexed
    return synthesize_timeseries(
        firings=F,
        waveforms=waveforms,
        timeseries_out=timeseries_out,
        noise_level=noise_level,
        samplerate=samplerate,
        duration=duration,
        waveform_upsamplefac=waveform_upsamplefac,
        amplitudes_row=amplitudes_row
    )
    
def time_basis_func(j,nnodes,times):
    if nnodes==1:
        return np.ones(times.shape)
    v=np.zeros(nnodes)
    v[j]=1
    tnodes=np.arange(nnodes)/(nnodes-1)
    if nnodes<3:
        kind0='linear'
    elif nnodes==3:
        kind0='quadratic'
    else:
        kind0='cubic'
    return spi.interp1d(tnodes,v,kind=kind0)(times)

def test_synthesize_drifting_timeseries():
    waveform_upsamplefac=5
    M=4 #num channels
    T=800 #num timepoints per clip
    K=2 # num units
    L=20
    num_interp_nodes=4
    from p_synthesize_random_waveforms import synthesize_random_waveforms
    #W=np.random.rand(M,T*waveform_upsamplefac,K*2)
    (W,geom)=synthesize_random_waveforms(M=M,T=T,K=K*num_interp_nodes,upsamplefac=waveform_upsamplefac)
    
    times=np.arange(L)/(L-1)*10000
    #labels=np.tile([1,2],int(L/2))
    labels=np.tile([1],L)
    F=np.array([np.zeros(L),times,labels])
    #return;
    X=synthesize_drifting_timeseries(waveforms=W,firings=F,duration=1,samplerate=10000,noise_level=0,num_interp_nodes=num_interp_nodes)
    import matplotlib.pyplot as pp
    pp.imshow(X,extent=[0, 1, 0, 1])
    #writemda32(X,'tmp.mda');
    return True

synthesize_drifting_timeseries.test=test_synthesize_drifting_timeseries
synthesize_drifting_timeseries.name = processor_name
synthesize_drifting_timeseries.version = processor_version

if __name__ == '__main__':
    print ('Running test')
    test_synthesize_drifting_timeseries()