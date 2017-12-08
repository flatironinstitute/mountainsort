import numpy as np

import sys, os

parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path+'/../mountainsort/packages/pyms')

from basic.p_compute_templates import compute_templates_helper
from basic.p_extract_clips import extract_clips_helper

from mlpy import readmda, writemda64, DiskReadMda

processor_name='pyms.apply_label_map'
processor_version='0.12'

def apply_label_map(*, firings, label_map, firings_out):
    """
    Apply a label map to a given firings, including masking and merging

    Parameters
    ----------
    firings : INPUT
        Path of input firings mda file
    label_map : INPUT
        Path of input label map mda file [base 1, mapping to zero removes from firings]
    firings_out : OUTPUT
        ...
    """
    firings =readmda(firings)
    label_map = readmda(label_map)
    label_map = np.reshape(label_map, (-1,2))
    label_map = label_map[np.argsort(label_map[:,0])] # Assure input is sorted

    #Propagate merge pairs to lowest label number
    for idx, label in enumerate(label_map[:,1]):
    	# jfm changed on 12/8/17 because isin() is not isin() older versions of numpy. :)
        #label_map[np.isin(label_map[:,0],label),0] = label_map[idx,0] # Input should be sorted
        label_map[np.where(label_map[:,0]==label)[0],0] = label_map[idx,0] # Input should be sorted

    #Apply label map
    for label_pair in range(label_map.shape[0]):
    	# jfm changed on 12/8/17 because isin() is not isin() older versions of numpy. :)
        #firings[2, np.isin(firings[2, :], label_map[label_pair, 1])] = label_map[label_pair,0]
        firings[2, np.where(firings[2, :] == label_map[label_pair, 1])[0]] = label_map[label_pair,0]

    #Mask out all labels mapped to zero
    firings = firings[:, firings[2, :] != 0]

    #Write remapped firings
    return writemda64(firings, firings_out)

apply_label_map.name = processor_name
apply_label_map.version = processor_version
apply_label_map.author = 'J Chung and J Magland'