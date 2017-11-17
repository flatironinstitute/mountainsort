import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from mlpy import writemda64,writemda32,readmda,DiskReadMda
from common import TimeseriesChunkReader

# we no longer use cppimport
# import cppimport
# cpp=cppimport.imp('basic_cpp')

# Do this first:
# g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` basic_cpp.cpp -o basic_cpp`python3-config --extension-suffix` -I../mlpy
import basic_cpp as cpp

processor_name='pyms.compute_templates'
processor_version='0.1'
def compute_templates(*,timeseries,firings,templates_out,clip_size=100):
    """
    Compute templates (average waveforms) for clusters defined by the labeled events in firings.

    Parameters
    ----------
    timeseries : INPUT
        Path of timeseries mda file (MxN) from which to draw the event clips (snippets) for computing the templates. M is number of channels, N is number of timepoints.
    firings : INPUT
        Path of firings mda file (RxL) where R>=3 and L is the number of events. Second row are timestamps, third row are integer labels.    
        
    templates_out : OUTPUT
        Path of output mda file (MxTxK). T=clip_size, K=maximum cluster label. Note that empty clusters will correspond to a template of all zeros. 
        
    clip_size : int
        (Optional) clip size, aka snippet size, number of timepoints in a single template
    """    
    templates=compute_templates_helper(timeseries=timeseries,firings=firings,clip_size=clip_size)
    return writemda32(templates,templates_out)
    
# Same as compute_templates, except return the templates as an array in memory
def compute_templates_helper(*,timeseries,firings,clip_size=100):
    X=DiskReadMda(timeseries)
    M,N = X.N1(),X.N2()
    N=N
    F=readmda(firings)
    L=F.shape[1]
    L=L
    T=clip_size
    times=F[1,:]
    labels=F[2,:].astype(int)
    K=np.max(labels)
    compute_templates._sums=np.zeros((M,T,K))
    compute_templates._counts=np.zeros(K)
    def _kernel(chunk,info):
        inds=np.where((info.t1<=times)&(times<=info.t2))[0]
        times0=(times[inds]-info.t1+info.t1a).astype(np.int32)
        labels0=labels[inds]
        
        clips0=np.zeros((M,clip_size,len(inds)),dtype=np.float32,order='F');
        cpp.extract_clips(clips0,chunk,times0,clip_size)
        
        for k in range(1,K+1):
            inds_kk=np.where(labels0==k)[0]
            compute_templates._sums[:,:,k-1]=compute_templates._sums[:,:,k-1]+np.sum(clips0[:,:,inds_kk],axis=2)
            compute_templates._counts[k-1]=compute_templates._counts[k-1]+len(inds_kk)
        return True
    TCR=TimeseriesChunkReader(chunk_size_mb=40, overlap_size=clip_size*2)
    if not TCR.run(timeseries,_kernel):
        return None
    templates=np.zeros((M,T,K))
    for k in range(1,K+1):
        if compute_templates._counts[k-1]:
            templates[:,:,k-1]=compute_templates._sums[:,:,k-1]/compute_templates._counts[k-1]
    return templates
    
compute_templates.name=processor_name
compute_templates.version=processor_version
def test_compute_templates():
    M,N,K,T,L = 5,1000,6,50,100
    X=np.random.rand(M,N)
    writemda32(X,'tmp.mda')
    F=np.zeros((3,L))
    F[1,:]=1+np.random.randint(N,size=(1,L))
    F[2,:]=1+np.random.randint(K,size=(1,L))
    writemda64(F,'tmp2.mda')
    ret=compute_templates(timeseries='tmp.mda',firings='tmp2.mda',templates_out='tmp3.mda',clip_size=T)
    assert(ret)
    templates0=readmda('tmp3.mda')
    assert(templates0.shape==(M,T,K))
    return True
compute_templates.test=test_compute_templates
