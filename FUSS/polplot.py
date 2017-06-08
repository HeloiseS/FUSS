"""
8 - June - 2017 / H.F.Stevance / fstevance1@sheffield.ac.uk

polplot.py is a submodule of the FUSS package. It is used to create polar plots.

Pre-requisites:
--------------
numpy, math, matplotlib.transforms, mpl_toolkits.axisartist.floating_axes, mpl_toolkits.axisartist.angle_helper,
matplotlib.projections, mpl_toolkits.axisartist.grid_finder

Functions:
----------
.axis(): Creates and returns a polar axis appropriate to the representation of SN ejecta. See docstring for me info.

.data(): Creates and returns the data to be plotted. Two 1D arrays of the same size are returned . The first contains
the angular component and the second contains the radial component of the data.

Working Example:
---------------

from FUSS import polplot as pp
import matplotlib.pyplot as plt

fig = plt.figure(1, figsize=(13, 6))

# List of locations and photospheric velocities for the 3 axes.
loc=[131, 132, 133]
vel=[20000, 18000, 15000]

# Plotting the 3 axes and storing them in the list g.
g=[]
for i in range(len(loc)):
    grid = pp.axis(fig, loc = loc[i], phot_vel=vel[i], ang_grid='l', rad_grid = 'ul', vel_lim=[0,40000])
    g.append(grid)

# Creating the data by giving the angle at start, angle at end, and velocity.
h = pp.data(15,30, 35000)
he = pp.data(20,40,22000)
ca = pp.data(120,130,33000)

# Plotting the data on the second axis (g[1])
g[1].plot(h[0], h[1], lw=3, label='H I')
g[1].plot(he[0], he[1], lw=3, label='He I')
g[1].plot(ca[0], ca[1], lw=3, label='Ca II')

# Creating a legend. Using bbox_to_anchor(x, y, width, height) so legend and plot do not overlap.
g[1].legend(bbox_to_anchor=(-3.5, 1.035, 5., .102), ncol=3)


plt.savefig('test_plot')
plt.show()

"""

import numpy as np
import math as m
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator, DictFormatter


def axis(fig, loc=111, num_ticks=1, phot_vel=None, vel_lim=[0, 30000], ang_grid ='ul', rad_grid = 'l'):
    """
    Sets up and returns the polar axis with no data. Data can be added later.
    :param fig: e.g fig = plt.figure(1, figsize=(10, 4))
    :param loc: Default = 111. Location of the subplot.
    :param num_ticks: Scales the number of ticks. Default is 1. Actual number of ticks is:
    int(maximum velocity/10000)*num_ticks
    :param phot_vel: Photospheric velocity. Give absolut value
    :param vel_lim: Velocity limits. Default [0, 30000]
    :param ang_grid: How many lines on angular grid. 'h' for heavy, 'l' for light, 'ul' for ultra-light.
    :param rad_grid: How many lines on radial grid. 'h' for heavy, 'l' for light, 'ul' for ultra-light, or specify where want the velocity grid using
    a list of velocities (again give absolute values, no negative numbers). Default is 'l'.
    :return: The polar plot! Now you can plot stuff on it.
    """

    # This will rotate the axes by 90 deg anti-clockwise
    tr_rotate = Affine2D().translate(90, 0)

    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    # scale degree to radians
    tr_scale = Affine2D().scale(np.pi/180., 1.)

    # This is the final translation, scaled to degrees, rotated and made polar.
    tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

    # This is the number of angular ticks, will change according to the type of grid chosen.
    if ang_grid == 'ul':
        grid_locator1 = angle_helper.LocatorD(3)
    if ang_grid == 'l':
        grid_locator1 = angle_helper.LocatorD(6)
    if ang_grid == 'h':
        grid_locator1 = angle_helper.LocatorD(10)

    # This is a tick formatter, I'm using it cuase the source code uses it and it works:
    tick_formatter1 = angle_helper.FormatterDMS()

    # Turns the radial tick labels from positive to negative. Maybe I could have just inverted radial axis
    # instead, but it was easier to copy those 2 lines from Emma's code since they work just fine.
    v_ticks=[(2500, '-2500'), (5000, '-5,000'), (7500, '-7,500'), (10000, '-10,000'), (12500, '-12,500'), (15000,'-15,000'), (17500, '-17,500'), (20000, '-20,000'),(22500, '-22500'), (25000, '-25,000'), (27500, '-27,500'), (30000, '-30,000'),(32500, '-32500'), (35000, '-35,000'), (37500, '-37,500'), (40000, '-40,000')] 
    tick_formatter2 = DictFormatter(dict(v_ticks))
    grid_locator2 = MaxNLocator(int(vel_lim[1]/10000)*num_ticks) # Set up number of ticks for the r-axis

    # the extremes are passed to the function
    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                extremes=(0, 180, vel_lim[0], vel_lim[1]),
                                grid_locator1=grid_locator1,
                                grid_locator2=grid_locator2,
                                tick_formatter1=tick_formatter1,
                                tick_formatter2=tick_formatter2,
                                )

    ax1 = floating_axes.FloatingSubplot(fig, loc, grid_helper=grid_helper)

    # Creating the subplot
    fig.add_subplot(ax1, pol=True)

    # The following are just to prettify the axes and their labels.
    # The axis artist lets you call axis with "bottom", "top", "left", "right". To know which one is
    # which just change it and try it to see what it does.

    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("top")

    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["left"].label.set_text("km/s")
    ax1.axis["left"].set_axis_direction("right")
    ax1.axis["top"].label.set_text(ur"$\theta$")

    # Creates a parasite axes whatever that is, but we need it.
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch  # Demo code comment: For aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9  # but this has a side effect that the patch is
                          # drawn twice, and possibly over some other
                          # artists. So, we decrease the zorder a bit to
                          # prevent this.
    
    # Angular markings, 'heavy' (every 15 deg), 'light' (every 30), or 'ultra-light' (every 45 deg).
    if ang_grid == 'h':
        mark_a=[15, 30, 45, 60, 75,  90, 105,  120, 135, 150, 165]
    if ang_grid == 'l':
        mark_a=[30, 60, 90, 120, 150]
    if ang_grid == 'ul':
        mark_a=[45,90,135]

    for a in mark_a:
        aux_ax.plot([a,a],[0,40000], color='0.35', linestyle='-.', linewidth=2)

    # Radial markings for the vel. Predefined "heavy", "light" or "ultra-light" option, or user can use their own list
    if rad_grid == 'h':
        mark_vel=[2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000, 27500,30000, 32500, 35000, 37500, 40000]  
    if rad_grid == 'l':
        mark_vel=[5000, 10000, 15000, 20000, 25000, 30000,35000,40000] 
    if rad_grid == 'ul':
        mark_vel = [10000, 20000, 30000,40000, 50000]     
    if rad_grid != 'h' and rad_grid != 'l' and rad_grid != 'ul':
        mark_vel = rad_grid
        
    for v in mark_vel:  
        y=np.arange(0, 180, 0.1)
        vaxis=np.empty(np.size(y))
        vaxis[:]=v
        aux_ax.plot(y,vaxis, color='0.35', linestyle=':', linewidth=1)

    if phot_vel != None:
        phot = datum(0, 180, phot_vel)
        aux_ax.plot(phot[0], phot[1], lw=3, c='k')
    
    return aux_ax


def datum(pa_start, pa_end, vel):
    """
    Creates one datum that can then be plotted on the polar axis or axes created.
    :param pa_start: Angle at start of range (in DEGREES)
    :param pa_end: Angle at end of range (in DEGREES)
    :param vel: Velocity
    :return: Two 1D arrays of the same length containing the angular component and radial component to be plotted
    """
    angle_range=np.arange(pa_start,pa_end+0.1,0.1)  #have to add 0.1 to pa_end to get full range

    # Making array of same length as angle_range and filling it wil the velocity given by user.
    vel1=np.empty(np.size(angle_range))
    vel1[:]=m.fabs(vel)

    return angle_range, vel1


def data(pa, pa_r, vel):
    pa_ranges=np.array([])
    vel_ranges=np.array([])

    for i in range(len(vel)):
        if pa_r[i] > 180:
            pa_r[i] = 180

        pa_start = pa[i] - pa_r[i]
        pa_end = pa[i] + pa_r[i]

        if pa_start < 0:
            print pa[i], pa_r[i]
            pa_start1 = 180-pa_start
            angle_range, vel1 = datum(pa_start1, 179.9, vel[i])
            pa_ranges = np.append(pa_ranges, [angle_range])
            vel_ranges = np.append(vel_ranges, [vel1])

            angle_range, vel1 = datum(0.1, pa_end, vel[i])
            pa_ranges = np.append(pa_ranges, [angle_range])
            vel_ranges = np.append(vel_ranges, [vel1])

        elif pa_end > 180:
            print pa[i], pa_r[i]
            pa_end1 = pa_end-180
            angle_range, vel1 = datum(pa_start, 179.9, vel[i])
            pa_ranges = np.append(pa_ranges, [angle_range])
            vel_ranges = np.append(vel_ranges, [vel1])

            angle_range, vel1 = datum(0.1, pa_end1, vel[i])
            pa_ranges = np.append(pa_ranges, [angle_range])
            vel_ranges = np.append(vel_ranges, [vel1])

        else:
            angle_range, vel1 = datum(pa_start, pa_end, vel[i])
            pa_ranges = np.append(pa_ranges, [angle_range])
            vel_ranges = np.append(vel_ranges, [vel1])

    return pa_ranges, vel_ranges