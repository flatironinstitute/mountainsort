import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import writemda32,writemda64,readmda,DiskReadMda,DiskWriteMda,MdaHeader
from common import TimeseriesChunkReader

processor_name='pyms.partial_timeseries'
processor_version='0.10'

def partial_timeseries(*,timeseries,timeseries_out_all=[],t1_all=[],t2_all=[],
                       timeseries_dtype='',timeseries_num_channels=4):
    """
    (Different from pyms.p_extract_timeseries, this function reads in a timeseries, a list of query data time points, 
        and output a list of cut files.)
    Extract a chunk of a timeseries dataset with events between A LIST OF t1's and t2's.

    Parameters
    ----------
    timeseries : INPUT
        Path of firing.mda timeseries (as it uses time ind in the second row to index how to cut), 
        If non-binary, data should be in shape of MxN where M is number of channels and N number of timepoints, in either .mda or raw binary format. 
        If raw binary, then you must supply dtype.
    timeseries_out_all: OUTPUT
        A list of paths of output timeseries in .mda format
    t1_all: list of integer, start timepoints (zero-based indexing).
    t2_all: list of integer, Integer end timepoints (zero-based indexing),
    timeseries_dtype : string
        Only supply this if timeseries is in raw binary format. Choices are int16, uint16, int32, float32, etc.
    timeseries_num_channels : integer
        Only supply this if timeseries is in raw binary format. Integer representing number of channels. Number of timepoints will be deduced
    """

    
    header0=None
    if (timeseries_dtype):
        size_bytes=os.path.getsize(timeseries)
        num_bytes_per_entry=get_num_bytes_per_entry_from_dt(timeseries_dtype)
        num_entries=size_bytes/num_bytes_per_entry
        if (num_entries % timeseries_num_channels != 0):
            print ("File size (%ld) is not divisible by number of channels (%g) for dtype=%s" % (size_bytes,timeseries_num_channels,timeseries_dtype))
            return False            
        num_timepoints=num_entries/timeseries_num_channels
        header0=MdaHeader(timeseries_dtype,[timeseries_num_channels,num_timepoints])
    
    X=DiskReadMda(timeseries,header0)
    M,N = X.N1(),X.N2()

    Xdata=X.readChunk(i1=0,i2=0,N1=M,N2=N)
        
    chunk_size_mb=100

    for i in range(0,len(t1_all)):
        t1=t1_all[i]
        t2=t2_all[i]
        timeseries_out=timeseries_out_all[i]

        t1t2=np.nonzero(np.logical_and([Xdata[1,:]>=t1],[Xdata[1,:]<=t2]))
        t1_update=t1t2[1][0]
        t2_update=t1t2[1][-1]
        N2=t2_update-t1_update+1 
        _writer=DiskWriteMda(timeseries_out,[M,N2],dt=X.dt())
        
        def _kernel(chunk,info):
            return _writer.writeChunk(chunk,i1=0,i2=info.t1)
        TCR=TimeseriesChunkReader(chunk_size_mb=chunk_size_mb, overlap_size=0, t1=t1_update, t2=t2_update)    
        TCR.run(X,_kernel)            
    
def get_num_bytes_per_entry_from_dt(dt):
	if dt == 'uint8':
		return 1
	if dt == 'float32':
		return 4
	if dt == 'int16':
		return 2
	if dt == 'int32':
		return 4
	if dt == 'uint16':
		return 2
	if dt == 'float64':
		return 8
	if dt == 'uint32':
		return 4
	return None

partial_timeseries.name=processor_name
partial_timeseries.version=processor_version  

def test_partial_timeseries():
    M,N = 4,12
    X=np.random.rand(M,N)
    X=np.ndarray((M,N))
    for n in range(0,N):
        for m in range(0,M):
            X[m,n]=n*10+m
    print('X',X)
    writemda32(X,'tmp1.mda')
    ret=partial_timeseries(timeseries="tmp1.mda",timeseries_out_all=["tmp_out1.mda","tmp_out2.mda"],t1_all=[20,21],t2_all=[80,82])
    #ret=extract_timeseries(timeseries="tmp.mda",timeseries_out="tmp2.mda",channels="1,3",t1=-1,t2=-1)

    A=readmda('tmp1.mda')
    B=readmda('tmp_out1.mda')
    C=readmda('tmp_out2.mda')
    print('A',A)
    print('tmp_out_from20_to80',B)
    print('tmp_out2_from21_to82',C)
    return True 
partial_timeseries.test=test_partial_timeseries

if __name__ == '__main__':
    print ('Running test')
    test_partial_timeseries()