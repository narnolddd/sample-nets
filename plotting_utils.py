import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def sqrt_custom(x1):
    if x1>0:
        return np.sqrt(x1)
    else:
        return 0.0

mysqrt = np.frompyfunc(sqrt_custom,1,1)

## shortcut to the plot/fill_inbetween fns
def plot_shaded_bars(xrange, means, sds, ax, marker, label=None):
    ax.plot(xrange,means,label=label, marker=marker)
    ax.fill_between(xrange,means - sds, means+sds, alpha=0.3)

def plot_shaded_bars_sqrt(xrange, means, sds, ax, marker, label=None):
    ax.plot(xrange,np.array(mysqrt(means)),label=label, marker=marker)
    ax.fill_between(xrange,np.array(mysqrt(means-sds),dtype=float), np.array(mysqrt(means+sds),dtype=float), alpha=0.3)