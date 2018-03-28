import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import readmda,writemda32

processor_name='pyms.synthesize_timeseries'
processor_version='0.11a'
def synthesize_timeseries(*,firings='',waveforms='',timeseries_out,noise_level=1,samplerate=30000,duration=60,waveform_upsamplefac=1,amplitudes_row=0):
    """
    Synthesize an electrophysiology timeseries from a set of ground-truth firing events and waveforms

    Parameters
    ----------
    firings : INPUT
        (Optional) The path of firing events file in .mda format. RxL where R>=3 and L is the number of events. Second row is the timestamps, third row is the integer labels/
    waveforms : INPUT
        (Optional) The path of (possibly upsampled) waveforms file in .mda format. Mx(T*waveform_upsample_factor)*K, where M is the number of channels, T is the clip size, and K is the number of units.
    
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
    """
    num_timepoints=np.int64(samplerate*duration)
    waveform_upsamplefac=int(waveform_upsamplefac)
    
    if type(waveforms)==str:
        if waveforms:
            W=readmda(waveforms)
        else:
            W=np.zeros((4,100*waveform_upsamplefac,0))
    else:
        W=waveforms
        
    
    
    if type(firings)==str:
        if firings:
            F=readmda(firings)
        else:
            F=np.zeros((3,0))
    else:
        F=firings
        
        
    times=F[1,:]
    labels=F[2,:].astype('int')
    
    M,TT,K = W.shape[0],W.shape[1],W.shape[2]
    T=int(TT/waveform_upsamplefac)
    Tmid=int(np.ceil((T+1)/2-1))
    
    N=num_timepoints
    if (N==0):
        if times.size==0:
            N=T
        else:
            N=max(times)+T
            
    X=np.random.randn(M,N)*noise_level

    waveform_list=[]
    for k in range(K):
        waveform0=W[:,:,k-1]
        waveform_list.append(waveform0)

    for j in range(times.size):
        t0=times[j]
        k0=labels[j]
        amp0=1
        if amplitudes_row>0:
            amp0=F[amplitudes_row-1,j]
        waveform0=waveform_list[k0-1]
        frac_offset=int(np.floor((t0-np.floor(t0))*waveform_upsamplefac))
        tstart=np.int64(np.floor(t0))-Tmid
        if (0<=tstart) and (tstart+T<=N):
            X[:,tstart:tstart+T]=X[:,tstart:tstart+T]+waveform0[:,frac_offset::waveform_upsamplefac]*amp0

    if timeseries_out:
        return writemda32(X,timeseries_out)
    else:
        return (X)

def test_synthesize_timeseries():
    #ret=synthesize_timeseries(timeseries_out='tmp.mda',duration=60,samplerate=30000)
    waveform_upsamplefac=13
    W=np.random.rand(4,100*waveform_upsamplefac,6)  
    F=np.zeros((3,0))
    X=synthesize_timeseries(waveforms=W,firings=F,duration=60,samplerate=30000,amplitudes_row=4)
    print(X.shape)
    #import matplotlib.pyplot as pp
    #pp.imshow(X,extent=[0, 5, 0, .3]);
    #writemda32(X,'tmp.mda');
    return True

synthesize_timeseries.test=test_synthesize_timeseries
synthesize_timeseries.name = processor_name
synthesize_timeseries.version = processor_version

if __name__ == '__main__':
    print ('Running test')
    test_synthesize_timeseries()