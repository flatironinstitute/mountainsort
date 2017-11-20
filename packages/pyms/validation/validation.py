import os
import sys

parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

from mlpy import ProcessorManager

import p_compute_accuracies

PM=ProcessorManager()

PM.registerProcessor(p_compute_accuracies.compute_accuracies)

if not PM.run(sys.argv):
    exit(-1)
