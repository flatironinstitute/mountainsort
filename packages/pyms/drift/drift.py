import os
import sys

parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import ProcessorManager

import p_concatenate_firings
import p_handle_drift_in_segment
import p_join_segments

PM=ProcessorManager()

PM.registerProcessor(p_concatenate_firings.concatenate_firings)
PM.registerProcessor(p_handle_drift_in_segment.handle_drift_in_segment)
PM.registerProcessor(p_join_segments.join_segments)

if not PM.run(sys.argv):
    exit(-1)
