#!/usr/bin/env python
import sys,os
import time
import pylab as plt
from sbmtm import sbmtm
import graph_tool.all as gt
import numpy as np
from matplotlib import pyplot as plt

gt.openmp_set_num_threads(int(sys.argv[1])) #set num threads
gt.seed_rng(42) #same results

print("Welcome to Topic Modelling")
print("using ",gt.openmp_get_num_threads(), " threads")

if __name__ == '__main__':
	start = time.time()
	print("initialised")
	gt.seed_rng(42)
	print("seed set")
	model = sbmtm()
	print("model created")
	model.load_graph(filename = 'graph.xml.gz')
	print("graph loaded")
	print(model.g)
	model.fit(n_init=2, verbose=True)
	model.savedata()
    print("mdl: %f"%model.get_mdl())
	print("ended")
	print("it took ", time.time()-start, "seconds")
