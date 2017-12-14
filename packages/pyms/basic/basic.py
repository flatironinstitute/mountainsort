import os
import sys

parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import ProcessorManager

import p_extract_geom
import p_extract_clips
import p_compute_templates
import p_extract_timeseries
#import p_bandpass_filter
import p_normalize_channels

PM=ProcessorManager()

PM.registerProcessor(p_extract_geom.extract_geom)
PM.registerProcessor(p_extract_clips.extract_clips)
PM.registerProcessor(p_compute_templates.compute_templates)
PM.registerProcessor(p_extract_timeseries.extract_timeseries)

# There is something wrong with bandpass filter. disabling for now
#PM.registerProcessor(p_bandpass_filter.bandpass_filter)


PM.registerProcessor(p_normalize_channels.normalize_channels)

if not PM.run(sys.argv):
    exit(-1)
