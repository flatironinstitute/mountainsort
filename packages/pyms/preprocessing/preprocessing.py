import sys

from mlpy import ProcessorManager

import p_bandpass_filter
import p_whiten

PM=ProcessorManager()

PM.registerProcessor(p_bandpass_filter.bandpass_filter)
PM.registerProcessor(p_whiten.whiten)

if not PM.run(sys.argv):
    exit(-1)
