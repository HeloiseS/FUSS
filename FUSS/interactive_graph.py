import matplotlib.pyplot as plt
import numpy as np
import FUSS as F


class DataRange():
    """
    I use this class to store data within a range.

    Notes
    ------
    Didn't put in a attribute for errors on x as for my spectra I usually don't have any.

    Attributes
    ----------
    name : string
        Name of the range
    x : int or float
        x coordinates of the data
    y : int or float
        y coordinates of the data
    yr : int or float
        errors on y
    start : int
        Beginning of the range (x dimension)
    end : int
        End of the range (x dimension)
    middle : int or float
        Median of the x coordinates
    avg: int or flaot
        Average of the y coordinates

    Methods
    -------
    average()

    """

    def __init__(self, name, xdata, ydata, ydata_err=None):
        """
        Examples
        --------
        instance_name = DataRange(name, xdata ydata, ydata_err=None)

        Notes
        -----
        xdata, ydata and ydata_err (if provided) must have the same dimensions

        Parameters
        ----------
        name : string
        xdata : array
            1D array containing the x coordinates of the data
        ydata : array
            1D array containing the y coordinates of the data
        ydata_err : array, optional
            1D array containing the error on the y coordinate of the data. Default is None.
            If errors are given the average function will perform a weighted average.

        """
        self.name = name
        self.x = xdata
        self.y = ydata
        self.yr = ydata_err
        self.start = min(self.x)
        self.end = max(self.x)
        self.middle = None
        self.avg = None

    def average(self):
        """
        Finds the average of ydata values and the median of xdata values and defines self.middle and self.avg.
        """
        if self.yr is None:
            self.middle = np.median(self.x)
            self.avg = np.average(self.y)
        else:
            self.middle = np.median(self.x)
            self.avg = np.average(self.y, weights=1 / (self.yr ** 2))


def onclick(event):
    """
    Takes in the  x coordinates of the point where the mouse has been clicked.
    Notes
    -----
    Requires ranges_graph (list), num (int) and coords (list), which are defined in def_ranges().

    Parameters
    ----------
    event : ?? e.g 'button_press_event'

    """
    global num
    global datx
    global coords  # stores x coordinates
    global ranges_graph

    datx = event.xdata
    coords.append(datx)
    if len(coords) == 2:
        # when 2 x values have been stored in coords they define a range
        # as such they are stored in ranges_graph. num is given too, this variable just numbers the ranges
        # that are in ranges_graph.
        print coords
        ranges_graph.append([num, coords[0], coords[1]])
        num += 1
        coords = []
    return


def def_ranges(fig, flux, err=False):
    """
    Defines ranges  of coordinates from user interaction with graph (mouse click).

    Notes
    -----
    Matplotlib figure should already be given and spectrum plotted created but not shown.

    Parameters
    ----------
    fig :
        matplotlib figure on which graph will be shown
    flux : array
        Flux spectrum data. 2-3 D array. flux[0] = wavelength, flux[1] = flux values. If err=True, flux[2] = error on
        flux[1]. Could be any 2-3 D data set.
    err : bool, optional
        Set to True if data has errors on y.

    Returns
    -------
    list of DataRange objects

    """
    # need to make some variables and lists globally defined for "onclick" to be able to use them
    global coords
    global ranges_graph
    global num
    coords = []
    ranges_graph = []
    num = 1

    cid = fig.canvas.mpl_connect('button_press_event', onclick)  # connects the button press event
    # on the figure fig to the onclick function.

    plt.show()  # shows graph, user shuld then click on graph to define x ranges

    print ranges_graph

    # the following takes the x ranges previously defined and create DataRange instances from the x, y and error
    # on y (if applicable) values for each range. These DataRange objects are then stored in ranges_data to be returned
    ranges_data = []
    for arange in ranges_graph:
        name = arange[0]  # the name of the range is just its number
        cond = (flux[0] > arange[1]) & (flux[0] < arange[2])
        xdata = flux[0][cond]
        ydata = flux[1][cond]
        #xdata = flux[0][np.argwhere(flux[0] > arange[1])[0]:np.argwhere(flux[0] < arange[2])[-1]]
        #ydata = flux[1][np.argwhere(flux[0] > arange[1])[0]:np.argwhere(flux[0] < arange[2])[-1]]
        if err is True:
            #ydata_err = flux[2][np.argwhere(flux[0] > arange[1])[0]:np.argwhere(flux[0] < arange[2])[-1]]
            ydata_err = flux[2][cond]
        else:
            ydata_err = None

        ranges_data.append(DataRange(name, xdata, ydata, ydata_err))

    return ranges_data
