import sys

from mlpy import ProcessorManager

import p_ms3alg

PM=ProcessorManager()

PM.registerProcessor(p_ms3alg.ms3alg)

if not PM.run(sys.argv):
    exit(-1)
