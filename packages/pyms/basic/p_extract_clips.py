import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import writemda64,writemda32,readmda,DiskReadMda
from common import TimeseriesChunkReader

# import the C++ code

# we no longer use cppimport
# import cppimport
# cpp=cppimport.imp('basic_cpp')

# Do this first:
# g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` basic_cpp.cpp -o basic_cpp`python3-config --extension-suffix` -I../mlpy
import basic_cpp as cpp

processor_name='pyms.extract_clips'
processor_version='0.1'
def extract_clips(*,timeseries,firings,clips_out,clip_size=100):
    """
    Extract clips corresponding to spike events

    Parameters
    ----------
    timeseries : INPUT
        Path of timeseries mda file (MxN) from which to draw the event clips (snippets)
    firings : INPUT
        Path of firings mda file (RxL) where R>=2 and L is the number of events. Second row are timestamps.
        
    clips_out : OUTPUT
        Path of clips mda file (MxTxL). T=clip_size
        
    clip_size : int
        (Optional) clip size, aka snippet size, aka number of timepoints in a single clip
    """    
    F=readmda(firings)
    times=F[1,:]
    clips=extract_clips_helper(timeseries=timeseries,times=times,clip_size=clip_size)
    return writemda32(clips,clips_out)

def extract_clips_helper(*,timeseries,times,clip_size=100,verbose=False):
    X=DiskReadMda(timeseries)
    M,N = X.N1(),X.N2()
    L=times.size
    T=clip_size
    extract_clips_helper._clips=np.zeros((M,T,L))
    def _kernel(chunk,info):
        inds=np.where((info.t1<=times)&(times<=info.t2))[0]
        times0=times[inds]-info.t1+info.t1a
        clips0=np.zeros((M,clip_size,len(inds)),dtype=np.float32,order='F');
        cpp.extract_clips(clips0,chunk,times0,clip_size)
        
        extract_clips_helper._clips[:,:,inds]=clips0
        return True
    TCR=TimeseriesChunkReader(chunk_size_mb=100, overlap_size=clip_size*2, verbose=verbose)
    if not TCR.run(timeseries,_kernel):
        return None
    return extract_clips_helper._clips
    

extract_clips.name=processor_name
extract_clips.version=processor_version
def test_extract_clips():
    M,T,L,N = 5,100,100,1000
    X=np.random.rand(M,N).astype(np.float32)
    writemda32(X,'tmp.mda')
    F=np.zeros((2,L))
    F[1,:]=200+np.random.randint(N-400,size=(1,L))
    writemda64(F,'tmp2.mda')
    ret=extract_clips(timeseries='tmp.mda',firings='tmp2.mda',clips_out='tmp3.mda',clip_size=T)
    assert(ret)
    clips0=readmda('tmp3.mda')
    assert(clips0.shape==(M,T,L))
    t0=int(F[1,10])
    a=int(np.floor((T+1)/2-1))
    np.array_equal(clips0[:,:,10],X[:,t0-a:t0-a+T])
    #np.testing.assert_almost_equal(clips0[:,:,10],X[:,t0-a:t0-a+T],decimal=4)
    return True
extract_clips.test=test_extract_clips

if __name__ == '__main__':
    print ('Running test')
    test_extract_clips()