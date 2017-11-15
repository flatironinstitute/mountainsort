import os
import sys

parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import ProcessorManager

import p_create_label_map
import p_apply_label_map

PM=ProcessorManager()

PM.registerProcessor(p_create_label_map.create_label_map)
PM.registerProcessor(p_apply_label_map.apply_label_map)

if not PM.run(sys.argv):
    exit(-1)
