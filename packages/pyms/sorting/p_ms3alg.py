from mlpy import mdaio
import numpy as np
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
            print ('Processed {} of {} chunks...'.format(self.num_completed_chunks.value,self.num_chunks))

def get_channel_neighborhood(m,Geom,adjacency_radius):
    M=Geom.shape[0]
    deltas=Geom-np.tile(Geom[m,:],(M,1))
    distsqrs=np.sum(Geom**2,axis=1)
    inds=np.where(distsqrs<=adjacency_radius**2)
    inds=np.sort(inds)
    return inds.ravel()

def detect_on_channel(data,detect_threshold,detect_interval,detect_sign):
    # Adjust the data to accommodate the detect_sign
    # After this adjustment, we only need to look for positive peaks
    if detect_sign<0:
        data=data*(-1)
    elif detect_sign==0:
        data=np.abs(data)
    elif detect_sign>0:
        pass
        
    #An event at timepoint t is flagged if the following two criteria are met:
    # 1. The value at t is greater than the detection threshold (detect_threshold)
    # 2. The value at t is greater than the value at any other timepoint within plus or minus <detect_interval> samples
    
    # First split the data into segments of size detect_interval (don't worry about timepoints left over, we assume we have padding)
    N=len(data)
    S2=int(np.floor(N/detect_interval))
    N2=S2*detect_interval
    data2=np.reshape(data[0:N2],(S2,detect_interval))
    
    # Find the maximum on each segment (these are the initial candidates)
    max_inds2=np.argmax(data2,axis=1)
    max_inds=max_inds2+detect_interval*np.arange(0,S2)
    max_vals=data[max_inds]
    
    # The following two tests compare the values of the candidates with the values of the neighbor candidates
    # If they are too close together, then discard the one that is smaller by setting its value to -1
    # Actually, this doesn't strictly satisfy the above criteria but it is close
    # TODO: fix the subtlety
    max_vals[np.where((max_inds[0:-1]>=max_inds[1:]-detect_interval) & (max_vals[0:-1]<max_vals[1:]))]=-1
    max_vals[1+np.array(np.where((max_inds[1:]<=max_inds[0:-1]+detect_interval) & (max_vals[1:]<=max_vals[0:-1])))]=-1
    
    # Finally we use only the candidates that satisfy the detect_threshold condition
    times=max_inds[np.where(max_vals>=detect_threshold)]
    return times

def extract_clips(data,times,clip_size):
    M=data.shape[0]
    T=clip_size
    L=len(times)
    Tmid = int(np.floor((T + 1) / 2) - 1);
    clips=np.zeros((M,T,L),dtype='float32')
    for j in range(L):
        t1=times[j]-Tmid
        t2=t1+clip_size
        clips[:,:,j]=data[:,t1:t2]
    return clips
    
def detect_and_extract_clips(num):
    opts=g_opts
    chunk_size=opts['chunk_size']
    clip_size=opts['clip_size']
    padding=opts['clip_size']*5
    timeseries_path=opts['timeseries']
    adjacency_radius=opts['adjacency_radius']
    detect_threshold=opts['detect_threshold']
    detect_interval=opts['detect_interval']
    detect_sign=opts['detect_sign']
    Geom=g_geom
    X=mdaio.DiskReadMda(timeseries_path)
    M=X.N1()
    N=X.N2()
    
    t1=int(num*opts['chunk_size']) # first timepoint of the chunk
    t2=int(np.minimum(X.N2(),(t1+chunk_size))) # last timepoint of chunk (+1)
    s1=int(np.maximum(0,t1-padding)) # first timepoint including the padding
    s2=int(np.minimum(X.N2(),t2+padding)) # last timepoint (+1) including the padding
    
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
    
    times_list=[]
    clips_list=[]
    for m in range(M):
        neighborhood=get_channel_neighborhood(m,Geom,adjacency_radius)
        padded_central_channel=padded_chunk[m,:].ravel()
        neighborhood_padded_chunk=padded_chunk[neighborhood,:]
        times_m=detect_on_channel(padded_central_channel,detect_threshold,detect_interval,detect_sign) # event times relative to the padded chunk
        times_m=times_m[np.where((times_m>=t1-s1)&(times_m<t2-s1))]
        clips_m=extract_clips(neighborhood_padded_chunk,times_m,clip_size)
        times_list.append(times_m)
        clips_list.append(clips_m)

    ###########################################################################################
    # Now we wait until we are ready to append to the output file
    # Note that we need to append in order, thus the shared_data object
    ###########################################################################################
    g_shared_info_detect_and_extract_clips.reportChunkCompleted(num) # Report that we have completed this chunk
    while True:
        if num == g_shared_info_detect_and_extract_clips.lastAppendedChunk()+1:
            break
        time.sleep(0.005) # so we don't saturate the CPU unnecessarily
    
    # Append the filtered chunk (excluding the padding) to the output file
    for m in range(M):
        times_path=opts['times_path_prefix']+'-{}.mda'.format(m)
        clips_path=opts['clips_path_prefix']+'-{}.mda'.format(m)
        mdaio.appendmda(times_list[m],times_path)
        mdaio.appendmda(clips_list[m],clips_path)
    
    # Report that we have appended so the next chunk can proceed
    g_shared_info_detect_and_extract_clips.reportChunkAppended(num)

    # Print status if it has been long enough
    if g_shared_info_detect_and_extract_clips.elapsedTime()>4:
        g_shared_info_detect_and_extract_clips.printStatus()
        g_shared_info_detect_and_extract_clips.resetTimer()
    
def ms3alg(*,
        timeseries,geom='',
        firings_out,
        adjacency_radius=-1,detect_interval=10,detect_threshold=3,detect_sign=0,
        clip_size=50,
        consolidate_clusters='true',consolidation_factor=0.9,
        merge_across_channels='true',fit_stage='true',
        chunk_size=30000*10,num_processes=os.cpu_count()):
    """
    MountainSort spike sorting (version 3) - largely consistent with mountainsortalg.ms3alg

    Parameters
    ----------
    timeseries : INPUT
        MxN raw timeseries array (M = #channels, N = #timepoints)
    geom : INPUT
        Optional geometry file (.csv format)
        
    firings_out : OUTPUT
        Firings array channels/times/labels (3xL, L = num. events)
        
    adjacency_radius : float
        Radius of local sorting neighborhood, corresponding to the geometry file (same units). 0 means each channel is sorted independently. -1 means all channels are included in every neighborhood.
    detect_interval : int
        The minimum number of timepoints between adjacent spikes detected on the same channel.
    detect_threshold : float
        Threshold for event detection, corresponding to the input file. So if the input file is normalized to have noise standard deviation 1 (e.g., whitened), then this is in units of std. deviations away from the mean.
    detect_sign : int
        Use 1, -1, or 0 to detect positive peaks, negative peaks, or both, respectively
    clip_size : int
        Size of extracted clips or snippets, used throughout
    """
    
    tempdir=os.environ.get('ML_PROCESSOR_TEMPDIR')
    if not tempdir:
        print ('Warning: environment variable ML_PROCESSOR_TEMPDIR not set. Using current directory.')
        tempdir='.'
    print ('Using tempdir={}'.format(tempdir))
    
    # Read the header of the timeseries input to get the num. channels and num. timepoints
    X=mdaio.DiskReadMda(timeseries)
    M=X.N1() # Number of channels
    N=X.N2() # Number of timepoints
    
    # Read the geometry file
    if geom:
        Geom = np.genfromtxt(geom, delimiter=',')
    else:
        Geom = np.zeros((M,2))
        
    if Geom.shape[0] != M:
        raise Exception('Incompatible dimensions between geom and timeseries: {} != {}'.format(Geom.shape[1],M))
    
    # For now, let's sort one neighborhood at a time
    #for neigh_ch in range(M):
    
    num_chunks=int(np.ceil(N/chunk_size))
    print ('Chunk size: {}, Num chunks: {}, Num processes: {}'.format(chunk_size,num_chunks,num_processes))
    
    opts={
        "timeseries":timeseries,
        "clip_size":clip_size,
        "adjacency_radius":adjacency_radius,
        "detect_threshold":detect_threshold,
        "detect_interval":detect_interval,
        "detect_sign":detect_sign,
        "chunk_size":chunk_size,
        "num_processes":num_processes,
        "num_chunks":num_chunks,
        "times_path_prefix":tempdir+'/times',
        "clips_path_prefix":tempdir+'/clips'
    }
    
    global g_opts
    g_opts=opts
    
    global g_geom
    g_geom=Geom
    
    global g_shared_info_detect_and_extract_clips
    g_shared_info_detect_and_extract_clips=SharedChunkInfo(num_chunks)
    
    for m in range(M):
        mdaio.writemda32(np.zeros([0]),opts['times_path_prefix']+'-{}.mda'.format(m))
        mdaio.writemda32(np.zeros([M,clip_size,0]),opts['clips_path_prefix']+'-{}.mda'.format(m))
    
    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(detect_and_extract_clips,range(num_chunks),chunksize=1)

    for m in range(M):
        clips_m=mdaio.readmda(opts['clips_path_prefix']+'-{}.mda'.format(m))
        print(clips_m.shape)
        
    return False
ms3alg.name='pyms.ms3alg-dont-use-yet'
ms3alg.version='0.1'