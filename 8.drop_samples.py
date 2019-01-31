#!python

# Drops samples from the VDS and verify

from hail import *
hc = HailContext(log = 'log/08.drop_samples.log', tmp_dir = 'tmp/hail')

vds = hc.read('../cohort.annotated.nHetHom.vds')

vds = vds.drop_samples()

vds.write('../cohort.annotated.nHetHom.no_samples.vds')

