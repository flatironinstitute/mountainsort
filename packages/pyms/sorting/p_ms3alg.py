from mlpy import mdaio
import numpy as np
import multiprocessing
import time
import os
import isosplit5
import sys

# import h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py
warnings.resetwarnings()

class SharedInfo():
    def __init__(self,num_channels):
        self.timer_timestamp = multiprocessing.Value('d',time.time(),lock=False)
        self.num_channels=num_channels
        self.num_completed_channels = multiprocessing.Value('l',0,lock=False)
        self.lock = multiprocessing.Lock()
    def acquireLock(self):
        self.lock.acquire()
    def releaseLock(self):
        self.lock.release()
    def reportChannelCompleted(self,num):
        self.num_completed_channels.value+=1
    def resetTimer(self):
        self.timer_timestamp.value=time.time()
    def elapsedTime(self):
        return time.time()-self.timer_timestamp.value
    def printStatus(self):
        print ('Processed {} of {} channels...'.format(self.num_completed_channels.value,self.num_channels))

def get_channel_neighborhood(m,Geom,adjacency_radius):
    M=Geom.shape[0]
    if adjacency_radius<0:
        return np.arange(M)
    deltas=Geom-np.tile(Geom[m,:],(M,1))
    distsqrs=np.sum(deltas**2,axis=1)
    inds=np.where(distsqrs<=adjacency_radius**2)[0]
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
    max_vals[ np.where((max_inds[0:-1]>=max_inds[1:]-detect_interval) & (max_vals[0:-1]<max_vals[1:]))[0] ]=-1
    max_vals[1+np.array( np.where((max_inds[1:]<=max_inds[0:-1]+detect_interval) & (max_vals[1:]<=max_vals[0:-1]))[0] )]=-1
    
    # Finally we use only the candidates that satisfy the detect_threshold condition
    times=max_inds[ np.where(max_vals>=detect_threshold)[0] ]
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
    
def prepare_timeseries_hdf5(timeseries_fname,timeseries_hdf5_fname,chunk_size,padding):
    chunk_size_with_padding=chunk_size+2*padding
    with h5py.File(timeseries_hdf5_fname,"w") as f:
        X=mdaio.DiskReadMda(timeseries_fname)
        M=X.N1() # Number of channels
        N=X.N2() # Number of timepoints
        num_chunks=int(np.ceil(N/chunk_size))
        f.create_dataset('chunk_size',data=[chunk_size])
        f.create_dataset('num_chunks',data=[num_chunks])
        f.create_dataset('padding',data=[padding])
        for j in range(num_chunks):
            padded_chunk=np.zeros((X.N1(),chunk_size_with_padding),dtype=X.dt())    
            t1=int(j*chunk_size) # first timepoint of the chunk
            t2=int(np.minimum(X.N2(),(t1+chunk_size))) # last timepoint of chunk (+1)
            s1=int(np.maximum(0,t1-padding)) # first timepoint including the padding
            s2=int(np.minimum(X.N2(),t2+padding)) # last timepoint (+1) including the padding
            
            # determine aa so that t1-s1+aa = padding
            # so, aa = padding-(t1-s1)
            aa = padding-(t1-s1)
            padded_chunk[:,aa:aa+s2-s1]=X.readChunk(i1=0,N1=X.N1(),i2=s1,N2=s2-s1) # Read the padded chunk

            for m in range(M):
                f.create_dataset('part-{}-{}'.format(m,j),data=padded_chunk[m,:].ravel())

def extract_neighborhood_chunk_from_timeseries_hdf5(timeseries_hdf5_fname,*,channels,chunk_num):
    parts=[]
    with h5py.File(timeseries_hdf5_fname,"r") as f:
        for m in channels:
            part=np.array(f.get('part-{}-{}'.format(m,chunk_num)))
            parts.append(part)
    num_timepoints=len(parts[0])
    ret=np.zeros((len(channels),num_timepoints))
    for ich in range(len(channels)):
        ret[ich,:]=parts[ich]
    return ret

def subsample_array(X,target_num):
    if X.size==0:
        return X
    inds=np.where(np.random.uniform(0,1,(X.size))<=target_num/X.size)[0]
    return X[inds]

def compute_principal_components_from_clips(clips,num_components):
    A=clips.reshape((clips.shape[0]*clips.shape[1],clips.shape[2]))
    u,s,vt=np.linalg.svd(A)
    u=u[:,:num_components]
    u=u.reshape((clips.shape[0],clips.shape[1],num_components))
    return u

def should_use_template(template0, *, consolidation_factor, detect_sign, central_channel):
    peak_location_tolerance=10
    M=template0.shape[0]
    T=template0.shape[1]
    Tmid=int((T+1)/2)-1

    if detect_sign<0:
        template0=template0*(-1)
    elif detect_sign==0:
        template0=np.abs(template0)
    elif detect_sign>0:
        pass

    peak_on_central_channel = 0;
    peak_t_on_central_channel = Tmid;
    peak_t_on_central_channel=np.argmax(template0[central_channel,:])
    peak_on_central_channel=template0[central_channel,peak_t_on_central_channel]
    peak_on_all_channels=np.max(template0)

    if peak_on_central_channel < peak_on_all_channels*consolidation_factor:
        return False
    if np.abs(peak_t_on_central_channel - Tmid) > peak_location_tolerance:
        return False
    return True

def consolidate_clusters_based_on_neighborhood_templates(templates,*,consolidation_factor,detect_sign,central_channel):
    M=templates.shape[0]
    T=templates.shape[1]
    K=templates.shape[2]
    last_label=0
    label_map=np.zeros((K+1,),dtype='int64')
    for k in range(1,K+1):
        template0=templates[:,:,k-1].reshape((M,T))
        if should_use_template(template0,consolidation_factor=consolidation_factor,detect_sign=detect_sign,central_channel=central_channel):
            label_map[k]=last_label+1
            last_label+=1
    return label_map

def apply_label_map(label_map,times,labels):
    new_labels=label_map[labels]
    inds_to_use=np.where(new_labels>0)[0]
    new_times=times[inds_to_use]
    new_labels=new_labels[inds_to_use]
    return (new_times,new_labels)

def compute_neighborhood_templates(times,labels,timeseries_hdf5_fname,*,neighborhood_channels,clip_size):
    K=int(np.max(labels))
    M_neigh=len(neighborhood_channels)
    clip_sums=np.zeros((M_neigh,clip_size,K),dtype='float64')
    clip_counts=np.zeros((K),dtype='float64')
    with h5py.File(timeseries_hdf5_fname,'r') as f:
        num_chunks=np.array(f.get('num_chunks'))[0]
        chunk_size=np.array(f.get('chunk_size'))[0]
        padding=np.array(f.get('padding'))[0]
    for chunk_num in range(num_chunks):
        padded_chunk=extract_neighborhood_chunk_from_timeseries_hdf5(timeseries_hdf5_fname,channels=neighborhood_channels,chunk_num=chunk_num)

        t1=chunk_num*chunk_size
        t2=(chunk_num+1)*chunk_size
        inds=np.where((t1<=times)&(times<t2))[0]
        times0=times[inds]-(t1-padding)
        labels0=labels[inds]

        # Extract the clips
        clips0=extract_clips(padded_chunk,times0,clip_size)
        for k in range(1,K+1):
            inds_k=np.where(labels0==k)[0]
            if len(inds_k)>0:
                clip_counts[k-1]+=len(inds_k)
                clip_sums[:,:,k-1]+=np.sum(clips0[:,:,inds_k],axis=2).reshape((M_neigh,clip_size))
    templates=np.zeros((M_neigh,clip_size,K))
    for k in range(1,K+1):
        if clip_counts[k-1]:
            templates[:,:,k-1]=clip_sums[:,:,k-1]/clip_counts[k-1]
    return templates

def branch_isosplit(features,*,depth=1):
    min_size_to_try_split=20
    labels1=isosplit5.isosplit5(features).ravel().astype('int64')
    if np.min(labels1)<0:
        raise Exception('Unexpected error in isosplit5.')
    K=int(np.max(labels1))
    if K<=1 or depth<=1:
        return labels1
    label_offset=0
    labels_new=np.zeros(labels1.shape,dtype='int64')
    for k in range(1,K+1):
        inds_k=np.where(labels1==k)[0]
        if len(inds_k)>min_size_to_try_split:
            labels_k=branch_isosplit(features[:,inds_k],depth=depth-1)
            K_k=int(np.max(labels_k))
            if K_k>10:
                mdaio.writemda32(features[:,inds_k],'example-{}.mda'.format(K_k))
            labels_new[inds_k]=label_offset+labels_k
            label_offset+=K_k
        else:
            labels_new[inds_k]=label_offset+1
            label_offset+=1
    return labels_new

def compute_fit_scores(labels,clips,templates):
    M=clips.shape[0]
    T=clips.shape[1]
    L=clips.shape[2]
    clip_fits=templates[:,:,labels-1]
    diffs=clips-clip_fits
    ss1=np.sum(diffs.reshape((M*T,L))**2,axis=0)
    ss2=np.sum(clips.reshape((M*T,L))**2,axis=0)
    amount_ss_reduced=ss2-ss1
    return amount_ss_reduced

def run_fit_stage(chunk_num):
    opts=g_opts
    shared_info=g_shared_info
    timeseries_hdf5_fname=opts['timeseries_hdf5_fname']
    chunk_size=opts['chunk_size']
    clip_size=opts['clip_size']
    padding=opts['padding']
    T=clip_size
    Tmid = int(np.floor((T + 1) / 2) - 1);
    Geom=g_geom
    M_global=Geom.shape[0]

    shared_info.acquireLock()
    with h5py.File(opts['times_labels_hdf5_fname'],"r") as f:
        channels=np.array(f.get('channels'))
        times=np.array(f.get('times'))
        labels=np.array(f.get('labels'))
        templates=np.array(f.get('templates'))
    shared_info.releaseLock()

    t1=chunk_num*chunk_size-padding
    t2=t1+chunk_size+2*padding
    inds_chunk=np.where((t1+clip_size<=times)&(times<t2-clip_size))[0]
    channels=channels[inds_chunk]
    times=times[inds_chunk]-t1
    labels=labels[inds_chunk]

    L=len(times)

    event_codes=np.zeros((L,)) # 0 means undetermined, 1 means using, -1 means not using

    all_chans=list(range(M_global))
    padded_chunk=extract_neighborhood_chunk_from_timeseries_hdf5(timeseries_hdf5_fname,channels=all_chans,chunk_num=chunk_num)
    clips=extract_clips(padded_chunk,times,clip_size=clip_size)
    scores_that_need_updating=np.ones((L,)) # 1 means needs updating
    scores=np.zeros((L,)) # 

    something_changed=True
    while something_changed:
        something_changed=False
        inds_0=np.where(event_codes==0)[0]
        if len(inds_0)==0: # Nothing left to determine, we are done
            break
        L_0=len(inds_0)
        times_0=times[inds_0]
        labels_0=labels[inds_0]
        scores_0=scores[inds_0]
        inds_1=np.where(scores_that_need_updating[inds_0])[0]
        scores_0[inds_1]=compute_fit_scores(labels_0[inds_1],clips[:,:,inds_0[inds_1]],templates)
        ii1=0
        ii2=0
        for ii in range(len(scores_0)):
            if scores_0[ii]<=0:
                event_codes[inds_0[ii]]=-1 # do not use
                something_changed=True
            else:
                tt=times_0[ii]
                while times_0[ii1]<tt-Tmid:
                    ii1+=1
                while (ii2<L_0) and (times_0[ii2]<tt-Tmid+T):
                    ii2+=1
                max_local_score=np.max(scores_0[ii1:ii2])
                if scores_0[ii]==max_local_score:
                    event_codes[inds_0[ii]]=1 # do use it
                    something_changed=True
                    for jj in range(ii1,ii2):
                        t_offset=tt-times_0[jj]
                        template0=templates[:,:,int(labels_0[jj]-1)]
                        if t_offset<0:
                            clips[:,0:T+t_offset,inds_0[jj]]-=template0[:,-t_offset:] # subtract off fit
                        else:
                            clips[:,t_offset:,inds_0[jj]]-=template0[:,:T-t_offset] # subtract off fit

    inds_to_use=np.where((event_codes==1)&(padding<=times)&(times<padding+chunk_size))[0]
    channels=channels[inds_to_use]
    times=times[inds_to_use]+t1
    labels=labels[inds_to_use]

    shared_info.acquireLock()
    with h5py.File(opts['times_labels_hdf5_fname'],"a") as f:
        f.create_dataset('times-{}'.format(chunk_num),data=times)
        f.create_dataset('labels-{}'.format(chunk_num),data=labels)
        f.create_dataset('channels-{}'.format(chunk_num),data=channels)
    shared_info.releaseLock()

    #scores=compute_fit_scores(labels,clips,templates)
    #inds_to_use=np.where(scores>0)[0]
    
    #clips=clips[:,:,inds_to_use]

    #inds0=np.where((padding<=times)&(times<padding+chunk_size))[0]
    #times=times[inds0]
    #labels=labels[inds0]
    #channels=channels[inds0]
    #times+=t1
    
    #
    #L_after=len(times)
    #print ('L before/after = {}/{}'.format(L_before,L_after))

def sort_on_neighborhood(channel_num):
    opts=g_opts
    timeseries_hdf5_fname=opts['timeseries_hdf5_fname']
    temp_hdf5_fname=opts['temp_hdf5_fname_prefix']+'-{}'.format(channel_num)
    num_chunks=opts['num_chunks']
    chunk_size=opts['chunk_size']
    clip_size=opts['clip_size']
    padding=opts['padding']
    timeseries_path=opts['timeseries']
    adjacency_radius=opts['adjacency_radius']
    detect_threshold=opts['detect_threshold']
    detect_interval=opts['detect_interval']
    detect_sign=opts['detect_sign']
    consolidation_factor=opts['consolidation_factor']
    num_features=20
    Geom=g_geom
    M_global=Geom.shape[0]

    neighborhood_channels=get_channel_neighborhood(channel_num,Geom,adjacency_radius)
    print ('Channel {}. Neighborhood:'.format(channel_num),neighborhood_channels)
    central_channel_ind=np.where(neighborhood_channels==channel_num)[0][0]

    target_num_clips_for_pca=200 # number of trial clips for getting the PCA components
    target_num_clips_for_pca_per_chunk=target_num_clips_for_pca/num_chunks

    # Detect and get clips for pca
    clips_for_pca_list=[]
    total_num_events=0
    total_num_timepoints=0
    with h5py.File(temp_hdf5_fname,'a') as tempf:    
        for chunk_num in range(num_chunks):
            padded_chunk=extract_neighborhood_chunk_from_timeseries_hdf5(timeseries_hdf5_fname,channels=neighborhood_channels,chunk_num=chunk_num)

            times0=detect_on_channel(padded_chunk[central_channel_ind,:].ravel(),detect_threshold,detect_interval,detect_sign)
            times0=times0[ np.where((times0>=padding)&(times0<padded_chunk.shape[1]-padding))[0] ]
            total_num_events+=len(times0)

            ## Extract a sample of the clips
            times0_for_pca=subsample_array(times0,target_num_clips_for_pca_per_chunk)
            clips0_for_pca=extract_clips(padded_chunk,times0_for_pca,clip_size)
            clips_for_pca_list.append(clips0_for_pca)

            ## Increment counters and save the times
            total_num_timepoints+=padded_chunk.shape[1]-2*padding
            tempf.create_dataset('times-{}'.format(chunk_num),data=times0)
            ############################

    print ('Detected {} events over {} timepoints on channel {}'.format(total_num_events,total_num_timepoints,channel_num))
    clips_for_pca=np.concatenate(clips_for_pca_list,axis=2)
    M_neigh=clips_for_pca.shape[0]
    mdaio.writemda32(clips_for_pca,'test_clips_for_pca-{}.mda'.format(channel_num))

    principal_components=compute_principal_components_from_clips(clips_for_pca,num_features)
    pca_vectors=principal_components.reshape((M_neigh*clip_size,num_features)).transpose()

    # Compute the event features
    with h5py.File(temp_hdf5_fname,'a') as tempf:    
        for chunk_num in range(num_chunks):
            padded_chunk=extract_neighborhood_chunk_from_timeseries_hdf5(timeseries_hdf5_fname,channels=neighborhood_channels,chunk_num=chunk_num)

            ## Read the times
            times0=np.array(tempf.get('times-{}'.format(chunk_num)))

            ## extract clips ####
            clips0=extract_clips(padded_chunk,times0,clip_size)
            
            ## Compute and save the event features
            features0=pca_vectors @ clips0.reshape((clips0.shape[0]*clips0.shape[1],clips0.shape[2]))

            ## Save the features
            tempf.create_dataset('features-{}'.format(chunk_num),data=features0)

    # Collect the features, and do the clustering
    times_list=[] # for concatenation
    features_list=[] # for concatenation
    with h5py.File(temp_hdf5_fname,'a') as tempf:
        for chunk_num in range(num_chunks):
            times0=np.array(tempf.get('times-{}'.format(chunk_num)))
            times0+=chunk_num*chunk_size-padding
            times_list.append(times0)
            features0=np.array(tempf.get('features-{}'.format(chunk_num)))
            features_list.append(features0)
    features=np.concatenate(features_list,axis=1)

    mdaio.writemda32(features,'test_features-{}.mda'.format(channel_num))

    labels=branch_isosplit(features,depth=3)
    times=np.concatenate(times_list)

    # First round of consolidation
    ## todo: subsample the times/labels for computing of templates
    neighborhood_templates=compute_neighborhood_templates(times,labels,timeseries_hdf5_fname,neighborhood_channels=neighborhood_channels,clip_size=clip_size)
    consolidation_factor1=0.6 #conservatively low
    label_map=consolidate_clusters_based_on_neighborhood_templates(neighborhood_templates,consolidation_factor=consolidation_factor1,detect_sign=detect_sign,central_channel=central_channel_ind)
    times,labels=apply_label_map(label_map,times,labels)

    mdaio.writemda32(neighborhood_templates,'test-neighborhood-templates-{}.mda'.format(channel_num))

    # Compute global templates
    ## todo: subsample the times/labels for computing of global templates
    global_templates=compute_neighborhood_templates(times,labels,timeseries_hdf5_fname,neighborhood_channels=list(range(M_global)),clip_size=clip_size)

    K=int(np.max(labels))
    neighborhood_sorting_result={
        "times":times,
        "labels":labels,
        "K":K,
        "global_templates":global_templates
    }
    return neighborhood_sorting_result

def ms3alg(*,
        timeseries,geom='',
        firings_out,firings_nofit_out='',
        adjacency_radius=-1,detect_interval=10,detect_threshold=3,detect_sign=0,
        clip_size=50,
        consolidation_factor=0.9,
        merge_across_channels='true',fit_stage='true',
        chunk_size=30000*10,num_workers=os.cpu_count()):
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
    firings_nofit_out : OUTPUT
        Firings array before fit stage
        
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
    fit_stage : string
        Whether or not to perform the fit stage (true or false)
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
    
    num_chunks=int(np.ceil(N/chunk_size))
    print ('Chunk size: {}, Num chunks: {}, Num workers: {}'.format(chunk_size,num_chunks,num_workers))

    opts={
        "timeseries":timeseries,
        "timeseries_hdf5_fname":tempdir+'/timeseries.hdf5',
        "times_labels_hdf5_fname":tempdir+'/times_labels.hdf5',
        "temp_hdf5_fname_prefix":tempdir+'/temp',
        "clip_size":clip_size,
        "padding":clip_size*10,
        "adjacency_radius":adjacency_radius,
        "detect_threshold":detect_threshold,
        "detect_interval":detect_interval,
        "detect_sign":detect_sign,
        "chunk_size":chunk_size,
        "num_workers":num_workers,
        "num_chunks":num_chunks,
        "consolidation_factor":consolidation_factor
    }
    
    print ('Preparing timeseries...')
    prepare_timeseries_hdf5(opts['timeseries'],opts['timeseries_hdf5_fname'],opts['chunk_size'],opts['padding'])
    
    global g_shared_info
    g_shared_info=SharedInfo(M)
    global g_opts
    g_opts=opts
    global g_geom
    g_geom=Geom
    
    # Sort on neighborhoods
    print ('Sorting on neighborhoods...')
    pool = multiprocessing.Pool(processes=num_workers)
    neighborhood_sortings=pool.map(sort_on_neighborhood,range(M),chunksize=1)

    # Collect results
    print ('Collecting results...')
    # Improve the following (don't load everything in memory)
    all_channels_list=[] # for concatenation
    all_times_list=[] # for concatenation
    all_labels_list=[] # for concatenation
    all_global_templates_list=[] # for concatenation
    all_global_template_central_channels_list=[] # central channels for the templates, for concat
    K_all=0
    for m in range(M):
        times_m=neighborhood_sortings[m]['times']
        labels_m=neighborhood_sortings[m]['labels']
        K_m=neighborhood_sortings[m]['K']
        global_templates_m=neighborhood_sortings[m]['global_templates']
        labels_m+=K_all
        K_all+=K_m
        all_channels_list.append(np.ones(times_m.shape)*m)
        all_times_list.append(times_m)
        all_labels_list.append(labels_m)
        all_global_templates_list.append(global_templates_m)
        all_global_template_central_channels_list.append(np.ones((K_m,))*m)
    all_times=np.concatenate(all_times_list)
    all_labels=np.concatenate(all_labels_list)
    all_channels=np.concatenate(all_channels_list)
    sorted_inds=np.argsort(all_times)
    all_channels=all_channels[sorted_inds]
    all_times=all_times[sorted_inds]
    all_labels=all_labels[sorted_inds]
    all_global_templates=np.concatenate(all_global_templates_list,axis=2)
    all_global_template_central_channels=np.concatenate(all_global_template_central_channels_list)
    
    mdaio.writemda32(all_global_templates,'test-all-templates.mda')

    if firings_nofit_out:
        # prepare firings nofit
        print ('Preparing firings nofit...')
        L=len(all_times)
        firings=np.zeros((3,L))
        firings[0,:]=all_channels+1 # 1-based indexing for the output file
        firings[1,:]=all_times
        firings[2,:]=all_labels

        # Write the nofit output
        print ('Writing firings_nofit file...')
        mdaio.writemda64(firings,firings_nofit_out)

    with h5py.File(opts['times_labels_hdf5_fname'],"a") as f:
        f.create_dataset('times',data=all_times)
        f.create_dataset('labels',data=all_labels)
        f.create_dataset('channels',data=all_channels)
        f.create_dataset('templates',data=all_global_templates)

    # Fit stage on chunks
    if fit_stage=='true':
        print ('Fit stage...')
        pool2 = multiprocessing.Pool(processes=num_workers)
        #pool2 = multiprocessing.Pool(processes=1)
        pool2.map(run_fit_stage,range(num_chunks),chunksize=1)

        channels_list=[]
        times_list=[]
        labels_list=[]
        with h5py.File(opts['times_labels_hdf5_fname'],"a") as f:
            for chunk_num in range(num_chunks):
                channels0=np.array(f.get('channels-{}'.format(chunk_num)))
                times0=np.array(f.get('times-{}'.format(chunk_num)))
                labels0=np.array(f.get('labels-{}'.format(chunk_num)))
                channels_list.append(channels0)
                times_list.append(times0)
                labels_list.append(labels0)
        all_channels=np.concatenate(channels_list)
        all_times=np.concatenate(times_list)
        all_labels=np.concatenate(labels_list)

    # prepare firings
    print ('Preparing firings...')
    L=len(all_times)
    firings=np.zeros((3,L))
    firings[0,:]=all_channels+1 # 1-based indexing for the output file
    firings[1,:]=all_times
    firings[2,:]=all_labels

    # Write the output
    print ('Writing firings file...')
    mdaio.writemda64(firings,firings_out)

    return True
ms3alg.name='pyms.ms3alg-dont-use-yet'
ms3alg.version='0.1'