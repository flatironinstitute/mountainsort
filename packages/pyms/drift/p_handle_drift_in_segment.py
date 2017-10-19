import numpy as np
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
import sys
import os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)
sys.path.append(parent_path+'/basic')

import mlpy
#import synthesis
import basic
from p_compute_templates import compute_templates_helper
from p_extract_clips import extract_clips_helper
import itertools as it

processor_name='pyms.handle_drift_in_segment'
processor_version='0.1'

def handle_drift_in_segment(*,timeseries,firings,firings_out):
    """
    Handle drift in segment.

    Parameters
    ----------
    timeseries : INPUT
        Path to preprocessed timeseries from which the events are extracted from (MxN)
    firings : INPUT
        Path of input firings mda file
    firings_out : OUTPUT
        Path of output drift-adjusted firings mda file
        ...
    """
    subcluster_size = 500 # Size of subclusters for comparison of merge candidate pairs
    bin_factor = 10 # subcluster_size / bin_factor = numbins for hist
    corr_comp_thresh = 0.95 # Minimum correlation in templates to consider as merge candidate
    clip_size=50
    n_pca_dim=10
            
    ## compute the templates
    templates=compute_templates_helper(timeseries=timeseries,firings=firings,clip_size=clip_size)
    templates=np.swapaxes(templates,0,1)
    templates=np.swapaxes(templates,2,0) #Makes templates of form Clust x Chan x Clipsize
    firings=mlpy.readmda(firings)
    print('templates',templates.shape)

    ## Determine the merge candidate pairs based on correlation
    subflat_templates=np.reshape(templates,(templates.shape[0],-1)) #flatten templates from templates from M x N x L (Clust x Chan x Clipsize) to (clust x flat)
    pairwise_idxs=np.array(list(it.chain.from_iterable(it.combinations(range(templates.shape[0]),2)))) #Generates 1D Array of all poss pairwise comparisons of clusters ([0 1 2] --> [0 1 0 2 1 2])
    pairwise_idxs=pairwise_idxs.reshape(-1,2) #Reshapes array, from above to readable [[0,1],[0,2],[1,2]]
    pairwise_corrcoef=np.zeros(pairwise_idxs.shape[0]) #Empty array for all pairs correlation measurements
    for row in range(pairwise_idxs.shape[0]): #Calculate the correlation coefficient for each pair of flattened templates
        pairwise_corrcoef[row]=np.corrcoef(subflat_templates[:,pairwise_idxs[row,0]],subflat_templates[:,pairwise_idxs[row,1]])[1,0]
    pairs_for_eval=np.array(pairwise_idxs[pairwise_corrcoef>=corr_comp_thresh]) #Threshold the correlation array, and use to index the pairwise comparison array

    ## Loop through the pairs for comparison

    for pair_to_test in range(pairs_for_eval.shape[0]):  # Iterate through pairs that are above correlation comparison threshold

        ## Extract out the times and labels corresponding to the pair
        firings_subset = firings[:, np.isin(firings[2, :], pairs_for_eval[pair_to_test, :] + 1)]  # Generate subfirings of only events from given pair, correct for base 0 vs. 1 difference
        test_labels = firings_subset[2, :]  # Labels from the pair of clusters
        test_eventtimes = firings_subset[1, :]  # Times from the pair of clusters
        sort_indices = np.argsort(test_eventtimes)  # there's no strict guarantee the firing times will be sorted, so adding a sort step for safety
        test_labels = test_labels[sort_indices]
        test_eventtimes = test_eventtimes[sort_indices]

        ## find the subcluster times and labels
        subcluster_event_indices = find_random_paired_events(test_eventtimes, test_labels, subcluster_size)
        subcluster_times = test_eventtimes[subcluster_event_indices]
        subcluster_labels = test_labels[subcluster_event_indices]

        ## Extract the clips for the subcluster
        subcluster_clips=extract_clips_helper(timeseries=timeseries,times=subcluster_times,clip_size=clip_size)

        ## Compute the centroids and project the clips onto the direction of the line connecting the two centroids

        # PCA to extract features of clips (number dim = n_pca_dim);
        subcluster_clips = np.reshape(subcluster_clips,(subcluster_clips.shape[0], -1))  # Flatten clips for PCA (expects 2d array)
        dimenReduc = PCA(n_components=n_pca_dim, whiten=True)
        clip_features = dimenReduc.fit_transform(subcluster_clips)

        # Use label data to separate clips into two groups, and adjust for base 0 vs base 1 difference
        A_indices = np.isin(subcluster_labels, pairs_for_eval[pair_to_test, 0] + 1)
        B_indices = np.isin(subcluster_labels, pairs_for_eval[pair_to_test, 1] + 1)
        clip_features_A = clip_features[A_indices, :]
        clip_features_B = clip_features[B_indices, :]

        # Calculate centroid
        centroidA = np.mean(clip_features_A, axis=0)
        centroidB = np.mean(clip_features_B, axis=0)

        # Project points onto line
        V = centroidA - centroidB
        V = np.tile(V, (clip_features.shape[0], 1))
        clip_1d_projs = np.einsum('ij,ij->i', clip_features, V)

    print('looped through all pairs')
        ##Histogram points and test cut point


def find_temporally_proximal_events(times,labels,subcluster_size):
    FirstEventInPair_idx=np.flatnonzero(np.diff(labels)) #Take diff (subtract neighboring value in array) of LABELS to determine when two neighboring events come from different clusters
    time_differences=times[FirstEventInPair_idx+1]-times[FirstEventInPair_idx] #Calculate time diff between each event pair where the events in the pair come from diff clust
    idx_for_eval=FirstEventInPair_idx[np.argsort(time_differences)[0:subcluster_size]] #Sort based on this time difference
    
    ##Eliminate all duplicate events. If events with labels A B A are used, and both pair AB and BA are valid, only count B once
    idx_for_eval = np.unique(np.append(idx_for_eval,idx_for_eval+1))

    return idx_for_eval

def find_random_paired_events(times, labels, subcluster_size):
    A_labl = np.isin(labels, 1)
    B_labl = np.isin(labels, 2)
    A_rands = np.random.randint(0, sum(A_labl), subcluster_size)
    B_rands = np.random.randint(0, sum(B_labl), subcluster_size)
    A_idxs = np.nonzero(A_labl)[0]
    A_idxs = A_idxs[A_rands]
    B_idxs = np.nonzero(B_labl)[0]
    B_idxs = B_idxs[B_rands]

    ##Eliminate all duplicate events. If events with labels A B A are used, and both pair AB and BA are valid, only count B once
    idx_for_eval = np.unique(np.append(A_idxs, B_idxs))

    return idx_for_eval

"""
def test_handle_drift_in_segment():
    duration=600
    samplerate=30000
    waveform_upsamplefac=13
    synthesis.synthesize_random_waveforms(waveforms_out='waveforms.mda',geometry_out='geom.csv',upsamplefac=waveform_upsamplefac)
    synthesis.synthesize_random_firings(firings_out='firings_true.mda',K=20,samplerate=samplerate,duration=duration)
    synthesis.synthesize_timeseries(firings='firings_true.mda',waveforms='waveforms.mda',timeseries_out='raw.mda',waveform_upsamplefac=waveform_upsamplefac,samplerate=samplerate,duration=duration)
    basic.bandpass_filter(timeseries='raw.mda',timeseries_out='filt.mda',freq_min=300,freq_max=6000,freq_wid=1000,samplerate=samplerate)
    basic.normalize_channels(timeseries='filt.mda',timeseries_out='test.pre.mda')
    handle_drift_in_segment(timeseries='test.pre.mda',firings='firings_true.mda',firings_out='test_drift.firings.mda')

    #test_clips=np.reshape(np.array([i for i in range(18)]),(2,3,3))

handle_drift_in_segment.test=test_handle_drift_in_segment
"""

handle_drift_in_segment.name = processor_name
handle_drift_in_segment.version = processor_version
handle_drift_in_segment.author = 'J Chung and J Magland'

if __name__ == '__main__':
    print('Running test')
    test_handle_drift_in_segment()
