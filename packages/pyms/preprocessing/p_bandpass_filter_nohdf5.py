from mlpy import mdaio
import numpy as np
from scipy import special
import multiprocessing
import time
import os

class SharedChunkInfo():
    def __init__(self,num_chunks):
        self.timer_timestamp = multiprocessing.Value('d',time.time(),lock=False)
        self.last_appended_chunk = multiprocessing.Value('l',-1,lock=False)
        self.num_chunks=num_chunks
        self.num_completed_chunks = multiprocessing.Value('l',0,lock=False)
        self.lock = multiprocessing.Lock()
    def reportChunkCompleted(self,num):
        with self.lock:
            self.num_completed_chunks.value+=1
    def reportChunkAppended(self,num):
        with self.lock:
            self.last_appended_chunk.value=num
    def lastAppendedChunk(self):
        with self.lock:
            return self.last_appended_chunk.value
    def resetTimer(self):
        with self.lock:
            self.timer_timestamp.value=time.time()
    def elapsedTime(self):
        with self.lock:
            return time.time()-self.timer_timestamp.value
    def printStatus(self):
        with self.lock:
            print('Processed {} of {} chunks...'.format(self.num_completed_chunks.value,self.num_chunks))

def create_filter_kernel(N,samplerate,freq_min,freq_max,freq_wid=1000):
    # Matches ahb's code /matlab/processors/ms_bandpass_filter.m
    # improved ahb, changing tanh to erf, correct -3dB pts  6/14/16
    
    T = N / samplerate # total time
    df = 1 / T # frequency grid
    relwid = 3.0; # relative bottom-end roll-off width param, kills low freqs by factor 1e-5.
    
    k_inds=np.arange(0,N)
    k_inds=np.where(k_inds<=(N+1)/2,k_inds,k_inds-N)
    
    fgrid=df*k_inds
    absf=np.abs(fgrid)
    
    val=np.ones(fgrid.shape)
    if freq_min!=0:
        val=val*(1 + special.erf(relwid * (absf - freq_min) / freq_min)) / 2
        val=np.where(np.abs(k_inds)<0.1,0,val) # kill DC part exactly
    if freq_max!=0:
        val=val*(1 - special.erf((absf - freq_max) / freq_wid)) / 2;
    val=np.sqrt(val) # note sqrt of filter func to apply to spectral intensity not ampl
    return val
    
def filter_chunk(num):
    #print('Filtering {}'.format(num))
    opts=g_opts
    #print('Filtering chunk {} of {}'.format(num,opts['num_chunks']))
    in_fname=opts['timeseries'] # The entire (large) input file
    out_fname=opts['timeseries_out'] # The entire (large) output file
    samplerate=opts['samplerate']
    freq_min=opts['freq_min']
    freq_max=opts['freq_max']
    freq_wid=opts['freq_wid']
    chunk_size=opts['chunk_size']
    padding=opts['padding']
    
    X=mdaio.DiskReadMda(in_fname)
    
    chunk_size_with_padding=chunk_size+2*padding
    padded_chunk=np.zeros((X.N1(),chunk_size_with_padding),dtype='float32')    
    
    t1=int(num*opts['chunk_size']) # first timepoint of the chunk
    t2=int(np.minimum(X.N2(),(t1+chunk_size))) # last timepoint of chunk (+1)
    s1=int(np.maximum(0,t1-padding)) # first timepoint including the padding
    s2=int(np.minimum(X.N2(),t2+padding)) # last timepoint (+1) including the padding
    
    # determine aa so that t1-s1+aa = padding
    # so, aa = padding-(t1-s1)
    aa = padding-(t1-s1)
    padded_chunk[:,aa:aa+s2-s1]=X.readChunk(i1=0,N1=X.N1(),i2=s1,N2=s2-s1) # Read the padded chunk
    
    # Do the actual filtering with a DFT with real input
    padded_chunk=np.fft.rfft(padded_chunk) 
    # Subtract off the mean of each channel unless we are doing only a low-pass filter
    if freq_min!=0:
        for m in range(padded_chunk.shape[0]):
            padded_chunk[m,:]=padded_chunk[m,:]-np.mean(padded_chunk[m,:])
    kernel=create_filter_kernel(chunk_size_with_padding,samplerate,freq_min,freq_max,freq_wid)
    kernel=kernel[0:padded_chunk.shape[1]] # because this is the DFT of real data
    padded_chunk=padded_chunk*np.tile(kernel,(padded_chunk.shape[0],1))
    padded_chunk=np.fft.irfft(padded_chunk)
    
    ###########################################################################################
    # Now we wait until we are ready to append to the output file
    # Note that we need to append in order, thus the shared_data object
    ###########################################################################################
    g_shared_data.reportChunkCompleted(num) # Report that we have completed this chunk
    while True:
        if num == g_shared_data.lastAppendedChunk()+1:
            break
        time.sleep(0.005) # so we don't saturate the CPU unnecessarily
    
    # Append the filtered chunk (excluding the padding) to the output file
    mdaio.appendmda(padded_chunk[:,padding:padding+(t2-t1)],out_fname)
    
    # Report that we have appended so the next chunk can proceed
    g_shared_data.reportChunkAppended(num)

    # Print status if it has been long enough
    if g_shared_data.elapsedTime()>4:
        g_shared_data.printStatus()
        g_shared_data.resetTimer()
    
def bandpass_filter(*,
        timeseries,timeseries_out,
        samplerate,freq_min,freq_max,freq_wid=1000,
        padding=3000,chunk_size=3000*10,num_processes=os.cpu_count()):
    """
    Apply a bandpass filter to a multi-channel timeseries

    Parameters
    ----------
    timeseries : INPUT
        MxN raw timeseries array (M = #channels, N = #timepoints)
        
    timeseries_out : OUTPUT
        Filtered output (MxN array)
        
    samplerate : float
        The sampling rate in Hz
    freq_min : float
        The lower endpoint of the frequency band (Hz)
    freq_max : float
        The upper endpoint of the frequency band (Hz)
    freq_wid : float
        The optional width of the roll-off (Hz)
    """
    X=mdaio.DiskReadMda(timeseries)
    M=X.N1() # Number of channels
    N=X.N2() # Number of timepoints
    
    num_chunks=int(np.ceil(N/chunk_size))
    print('Chunk size: {}, Padding: {}, Num chunks: {}, Num processes: {}'.format(chunk_size,padding,num_chunks,num_processes))
    
    opts={
        "timeseries":timeseries,
        "timeseries_out":timeseries_out,
        "samplerate":samplerate,
        "freq_min":freq_min,
        "freq_max":freq_max,
        "freq_wid":freq_wid,
        "chunk_size":chunk_size,
        "padding":padding,
        "num_processes":num_processes,
        "num_chunks":num_chunks
    }
    global g_shared_data
    g_shared_data=SharedChunkInfo(num_chunks)
    global g_opts
    g_opts=opts
    mdaio.writemda32(np.zeros([M,0]),timeseries_out)
    
    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(filter_chunk,range(num_chunks),chunksize=1)
    return True
bandpass_filter.name='pyms.bandpass_filter'
bandpass_filter.version='0.1'