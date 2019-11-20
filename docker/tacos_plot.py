import numpy as np
from scipy.interpolate import interpn
import pandas as pd
import matplotlib.pyplot as plt

def scatterdense(x, y, ax=None, nbins=80, **kwargs):
    xmin = np.log10(1e-3)
    xmax = np.log10(x.max())
    ymin = np.log10(1e-6)
    ymax = np.log10(y.max())

    xbins = np.logspace(xmin, xmax, nbins) # <- make a range from 10**xmin to 10**xmax
    ybins = np.logspace(ymin, ymax, nbins) # <- make a range from 10**ymin to 10**ymax
    data , x_e, y_e = np.histogram2d(x, y, bins = (xbins, ybins))
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    if ax is None:
        fig=plt.figure()
        ax=fig.subplots()
    ax.scatter(x, y, c=z, **kwargs)