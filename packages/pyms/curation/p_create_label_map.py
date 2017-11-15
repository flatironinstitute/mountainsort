import numpy as np
import json
import sys, os

parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path+'/../mountainsort/packages/pyms')

from basic.p_compute_templates import compute_templates_helper
from basic.p_extract_clips import extract_clips_helper

from mlpy import readmda, writemda64, DiskReadMda

processor_name='pyms.create_label_map'
processor_version='0.1'

def create_label_map(*, metrics, label_map_out):
    """
    Generate a label map based on the metrics file, where labels being mapped to zero are to be removed.

    Parameters
    ----------
    metrics : INPUT
        Path of metrics json file to be used for generating the label map
    label_map_out : OUTPUT
        Path to mda file where the second column is the present label, and the first column is the new label
        ...
    """
    #TODO: Way to pass in logic or thresholds flexibly

    firing_rate_thresh = .05
    isolation_thresh = .95
    noise_overlap_thresh = .03
    peak_snr_thresh = 1.5

    label_map = []

    #Load json
    with open(metrics) as metrics_json:
        metrics_data = json.load(metrics_json)

    #Iterate through all clusters
    for idx in range(len(metrics_data['clusters'])):
        if metrics_data['clusters'][idx]['metrics']['firing_rate'] < firing_rate_thresh or \
            metrics_data['clusters'][idx]['metrics']['isolation'] < isolation_thresh or \
            metrics_data['clusters'][idx]['metrics']['noise_overlap'] > noise_overlap_thresh or \
            metrics_data['clusters'][idx]['metrics']['peak_snr'] < peak_snr_thresh:
            #Map to zero (mask out)
            label_map.append([0,metrics_data['clusters'][idx]['label']])
        else: #Check if burst parent exists
            if metrics_data['clusters'][idx]['metrics']['bursting_parent']:
                label_map.append([metrics_data['clusters'][idx]['metrics']['bursting_parent'],
                                  metrics_data['clusters'][idx]['label']])

    #Writeout
    return writemda64(np.array(label_map),label_map_out)

create_label_map.name = processor_name
create_label_map.version = processor_version
create_label_map.author = 'J Chung and J Magland'