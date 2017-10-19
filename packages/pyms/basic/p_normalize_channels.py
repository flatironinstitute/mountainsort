import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import writemda32,readmda,DiskReadMda,DiskWriteMda
from common import TimeseriesChunkReader

processor_name='pyms.normalize_channels'
processor_version='0.1'
def normalize_channels(*,timeseries,timeseries_out):
    """
    Normalize the channels in a timeseries array to each have unit variance

    Parameters
    ----------
    timeseries : INPUT
        Path of timeseries, MxN where M is number of channels and N number of timepoints, in .mda format
        
    timeseries_out : OUTPUT
        Path of output timeseries in .mda format            
    """
    
    X=DiskReadMda(timeseries)
    M,N = X.N1(),X.N2()    
    _writer=DiskWriteMda(timeseries_out,[M,N],dt=X.dt())
    
    chunk_size_mb=100
    normalize_channels._sums=np.zeros(M)
    normalize_channels._sumsqrs=np.zeros(M)

    def _kernel_compute_sumsqrs(chunk,info):
        normalize_channels._sums=normalize_channels._sums+np.sum(chunk,axis=1)
        normalize_channels._sumsqrs=normalize_channels._sumsqrs+np.sum(chunk**2,axis=1)
        return True
        
    def _kernel_normalize_and_write(chunk,info):
        Nchunk=chunk.shape[1]
        means=normalize_channels._sums/N
        variances=(normalize_channels._sumsqrs-normalize_channels._sums**2/N)/(N-1)
        stdevs=np.sqrt(variances)
        stdevs[np.where(stdevs==0)]=1
        means=np.reshape(means,(M,1))
        stdevs=np.reshape(stdevs,(M,1))
        chunk=(chunk-np.tile(means,(1,Nchunk)))/np.tile(stdevs,(1,Nchunk))
        return _writer.writeChunk(chunk,i1=0,i2=info.t1)
    
    TCR=TimeseriesChunkReader(chunk_size_mb=chunk_size_mb, overlap_size=0)    
    if not TCR.run(timeseries,_kernel_compute_sumsqrs):
        return False
    if not TCR.run(timeseries,_kernel_normalize_and_write):
        return False
    return True

normalize_channels.name=processor_name
normalize_channels.version=processor_version  
def test_normalize_channels():
    M,N = 4,1000
    X=np.random.rand(M,N)
    writemda32(X,'tmp.mda')
    ret=normalize_channels(timeseries="tmp.mda",timeseries_out="tmp2.mda")
    assert(ret)
    A=readmda('tmp.mda')
    B=readmda('tmp2.mda')
    A_mean=np.mean(A,axis=1)
    A_stdev=np.sqrt(np.var(A,axis=1,ddof=1))
    A_norm=(A-np.tile(np.reshape(A_mean,(M,1)),(1,N)))/np.tile(np.reshape(A_stdev,(M,1)),(1,N))
    np.testing.assert_array_almost_equal(A_norm,B,decimal=5)
    return True 
normalize_channels.test=test_normalize_channels

if __name__ == '__main__':
    print ('Running test')
    test_normalize_channels()