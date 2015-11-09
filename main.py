#!/usr/bin/env python

from simuseq import SimuSEQ

SimuSEQ().auto([1,1,1,1,0,1],'rear',10)
#setting[0] = SimuSEQ().readseq()
#setting[1] = SimuSEQ().makeprobe()
#setting[2] = SimuSEQ().microarray()
#setting[3] = SimuSEQ().rank()
#setting[4] = SimuSEQ().fasta()
#setting[5] = SimuSEQ().makedata()
#([setting],'front',probe_size) = front_probe
#([setting],'rear',probe_size) = rear_probe
#([setting],'full','full') = fullsize_probe
