import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)
sys.path.append(parent_path+'/basic')

#from p_compute_templates import compute_templates_helper
from p_extract_clips import extract_clips_helper

from mlpy import readmda,writemda64,DiskReadMda

processor_name='pyms.join_segments'
processor_version='0.1'
def join_segments(*,timeseries_list, firings_list, dmatrix_out, templates_out):
    """
    Join the results of spike sorting on a sequence of time segments to form a single firings file

    Parameters
    ----------
    timeseries_list : INPUT
        A list of paths of adjacent preprocessed timeseries segment files
    firings_list : INPUT
        A list of paths to corresponding firings files
        
    dmatrix_out : OUTPUT
        dmatrix for debugging    
    templates_out : OUTPUT
        templates for debugging

    """
    X=DiskReadMda(timeseries_list[0])
    M=X.N1()
    clip_size=100
    num_segments=len(timeseries_list)
    firings_arrays=[]
    for j in range(num_segments):
        F=readmda(firings_list[j])
        firings_arrays.append(F)
    Kmax=0;
    for j in range(num_segments):
        F=firings_arrays[j]
        labels=F[2,:]
        Kmax=int(max(Kmax,np.max(labels)))    
    dmatrix=np.ones((Kmax,Kmax,num_segments-1))*(-1)
    templates=np.zeros((M,clip_size,Kmax,2*(num_segments-1)))
    
    for j in range(num_segments-1):
        print ('Computing dmatrix between segments %d and %d' % (j,j+1))
        (dmatrix0,templates1,templates2)=compute_dmatrix(timeseries_list[j],timeseries_list[j+1],firings_arrays[j],firings_arrays[j+1],clip_size=clip_size)
        dmatrix[0:dmatrix0.shape[0],0:dmatrix0.shape[1],j]=dmatrix0
        templates[:,:,0:dmatrix0.shape[0],j*2]=templates1
        templates[:,:,0:dmatrix0.shape[1],j*2+1]=templates2
    
    writemda64(templates,templates_out)
    return writemda64(dmatrix,dmatrix_out)
    
def compute_dmatrix(timeseries1,timeseries2,F1,F2,*,clip_size):
    X=DiskReadMda(timeseries1)
    M=X.N1()
    F1b=get_last_events(F1,100)
    F2b=get_first_events(F2,100)
    times1=F1b[1,:].ravel()
    labels1=F1b[2,:].ravel()
    clips1=extract_clips_helper(timeseries=timeseries1,times=times1,clip_size=clip_size)
    times2=F2b[1,:].ravel()
    labels2=F2b[2,:].ravel()
    clips2=extract_clips_helper(timeseries=timeseries2,times=times2,clip_size=clip_size)

    K1=int(max(labels1))
    K2=int(max(labels2))
    dmatrix=np.zeros((K1,K2))
    templates1=np.zeros((M,clip_size,K1))
    templates2=np.zeros((M,clip_size,K2))
    for k1 in range(1,K1+1):
        #times1_k1=times1[np.where(labels1==k1)[0]]
        inds_k1=np.where(labels1==k1)[0]
        clips1_k1=clips1[:,:,inds_k1]
        templates1[:,:,k1-1]=np.mean(clips1_k1,axis=2)
        for k2 in range(1,K2+1):            
            #times2_k2=times2[np.where(labels2==k2)[0]]            
            inds_k2=np.where(labels2==k2)[0]            
            clips2_k2=clips2[:,:,inds_k2]            
            templates2[:,:,k2-1]=np.mean(clips2_k2,axis=2)
            dmatrix[k1-1,k2-1]=compute_distance_between_clusters(clips1_k1,clips2_k2)
    return (dmatrix,templates1,templates2)        
    
def get_first_events(firings,num):
    L=firings.shape[1]
    times=firings[1,:]
    labels=firings[2,:]
    K=int(max(labels))
    to_use=np.zeros(L)
    for k in range(1,K+1):
        inds_k=np.where(labels==k)[0]
        times_k=times[inds_k]
        if (len(times_k)<=num):
            to_use[inds_k]=1
        else:
            times_k_sorted=np.sort(times_k)
            cutoff=times_k_sorted[num]
            to_use[inds_k[np.where(times_k<=cutoff)[0]]]=1
    return firings[:,np.where(to_use==1)[0]]

def get_last_events(firings,num):
    L=firings.shape[1]
    times=firings[1,:]
    labels=firings[2,:]
    K=int(max(labels))
    to_use=np.zeros(L)
    for k in range(1,K+1):
        inds_k=np.where(labels==k)[0]
        times_k=times[inds_k]
        if (len(times_k)<=num):
            to_use[inds_k]=1
        else:
            times_k_sorted=np.sort(times_k)
            cutoff=times_k_sorted[len(times_k_sorted)-num]
            to_use[inds_k[np.where(times_k>=cutoff)[0]]]=1
    return firings[:,np.where(to_use==1)[0]]

    
def compute_distance_between_clusters(clips1,clips2):
    centroid1=np.mean(clips1,axis=2)
    centroid2=np.mean(clips2,axis=2)    
    dist=np.sum((centroid2-centroid1)**2)
    return dist

def test_join_segments():
    timeseries_list=['test1/pre_seg1.mda','test1/pre_seg2.mda','test1/pre_seg3.mda','test1/pre_seg4.mda','test1/pre_seg5.mda','test1/pre_seg6.mda']
    firings_list=['test1/firings_seg1.mda','test1/firings_seg2.mda','test1/firings_seg3.mda','test1/firings_seg4.mda','test1/firings_seg5.mda','test1/firings_seg6.mda']
    dmatrix_out='test1/dmatrix.mda'
    templates_out='test1/templates.mda'
    ret=join_segments(timeseries_list=timeseries_list,firings_list=firings_list,dmatrix_out=dmatrix_out,templates_out=templates_out)
    return ret

#join_segments.test=test_join_segments
join_segments.name = processor_name
join_segments.version = processor_version
join_segments.author = 'J Magland'

if __name__ == '__main__':
    print ('Running test')
    test_join_segments()