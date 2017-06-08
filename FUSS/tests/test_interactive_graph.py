import matplotlib.pyplot as plt
import FUSS.interactive_graph as ig
import numpy as np

def test_def_ranges():
    xdata = np.arange(0,10,0.1) # creates a dummy data set
    ydata = xdata*2
    fig = plt.figure() # need to plot a figure to give as input in ig.def_ranges
    plt.scatter(xdata, ydata) # plto the data here, plt.show() is called in def_ranges
    ranges = ig.def_ranges(fig, [xdata, ydata]) # use the function and retrieve the ranges 
                                                # this requires the user to click on plot when 
                                                # running the test
    assert len(ranges[0].x) < len(xdata) # here we check that the range returned is smaller 
                                         # than the whole data set. Should be the case if 
                                         # the user didn't select the whole data set. 
