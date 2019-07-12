#!/usr/bin/env python
import os
import pylab as plt
from sbmtm import sbmtm
import graph_tool.all as gt
import numpy as np
from matplotlib import pyplot as plt

print("welcome to Topic Modelling")

if __name__ == '__main__':
        print("initialised")
        gt.seed_rng(42)
        print("seed set")
        model = sbmtm()
        print("model created")
        model.load_graph(filename = 'graph.xml.gz')
        print("graph loaded")
        model.fit(verbose=True)
        model.savedata()
        print("ended")
