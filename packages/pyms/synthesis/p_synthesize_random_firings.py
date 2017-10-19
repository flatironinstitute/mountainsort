import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import readmda,writemda64

processor_name='pyms.synthesize_random_firings'
processor_version='0.13b'
def synthesize_random_firings(*,firings_out,K=20,samplerate=30000,duration=60):
    """
    Synthesize random waveforms for use in creating a synthetic timeseries dataset

    Parameters
    ----------
    firings_out : OUTPUT
        Path to output firings mda file. 3xL, L is the number of events, second row are timestamps, third row are integer unit labels
    
    K : int
        (Optional) number of simulated units
    samplerate : double
        (Optional) sampling frequency in Hz
    duration : double
        (Optional) duration of the simulated acquisition in seconds
    """
    firing_rates=3*np.ones((K))
    refr=10
    
    N=np.int64(duration*samplerate)
    
    # events/sec * sec/timepoint * N
    populations=np.ceil(firing_rates/samplerate*N).astype('int')
    times=np.zeros(0)
    labels=np.zeros(0)
    for k in range(1,K+1):
        refr_timepoints=refr/1000*samplerate
        
        times0=np.random.rand(populations[k-1])*(N-1)+1

        ## make an interesting autocorrelogram shape
        times0=np.hstack((times0,times0+rand_distr2(refr_timepoints,refr_timepoints*20,times0.size)))
        times0=times0[np.random.choice(times0.size,int(times0.size/2))]
        times0=times0[np.where((0<=times0)&(times0<N))]
        
        times0=enforce_refractory_period(times0,refr_timepoints)
        times=np.hstack((times,times0))
        labels=np.hstack((labels,k*np.ones(times0.shape)))

    sort_inds=np.argsort(times)
    times=times[sort_inds]
    labels=labels[sort_inds]

    firings=np.zeros((3,times.size),dtype=np.float64)
    firings[1,:]=times
    firings[2,:]=labels
    return writemda64(firings,firings_out)
    
def rand_distr2(a,b,num):
    X=np.random.rand(num)
    X=a+(b-a)*X**2
    return X
    
def enforce_refractory_period(times_in,refr):
    if (times_in.size==0): return times_in
    
    times0=np.sort(times_in)
    done=False
    while not done:
        diffs=times0[1:]-times0[:-1]
        diffs=np.hstack((diffs,np.inf)) #hack to make sure we handle the last one
        inds0=np.where((diffs[:-1]<=refr)&(diffs[1:]>=refr))[0] #only first violator in every group
        if len(inds0)>0:
            times0[inds0]=-1 #kind of a hack, what's the better way?
            times0=times0[np.where(times0>=0)]
        else:
            done=True
        
    return times0
        
def test_synthesize_random_firings():
    K=10
    synthesize_random_firings(K=K,firings_out='tmp.firings.mda')
    firings=readmda('tmp.firings.mda')
    labels=firings[2,:]
    assert(max(labels)==K)
    assert(firings.shape[0]==3)    
    return True

synthesize_random_firings.test=test_synthesize_random_firings
synthesize_random_firings.name = processor_name
synthesize_random_firings.version = processor_version

if __name__ == '__main__':
    print ('Running test')
    test_synthesize_random_firings()
    
