#!/usr/local/bin/python

from simuseq import SimuSEQ

seq = SimuSEQ().readseq()
SimuSEQ().makeprobe(seq)
SimuSEQ().microarray(seq)
SimuSEQ().fasta(seq)
SimuSEQ().makedata()
