import sys

from mlpy import ProcessorManager

import p_bandpass_filter

PM=ProcessorManager()

PM.registerProcessor(p_bandpass_filter.bandpass_filter)

if not PM.run(sys.argv):
    exit(-1)
