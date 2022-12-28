import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

def reset_plt(ticksize,fontsize):
    plt.style.use('seaborn-white')
    plt.rcParams['xtick.labelsize'] = ticksize
    plt.rcParams['ytick.labelsize'] = ticksize
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['legend.facecolor'] = 'white'
    plt.rcParams['axes.formatter.limits'] = (-1, 3)
    plt.rcParams['axes.linewidth'] = 2.25
    

def put_ticks(this_fig,this_ax):
    this_ax.xaxis.set_tick_params(which='major', direction='in', width=2.5, length=12, zorder=1, top=True)
    this_ax.yaxis.set_tick_params(which='major', direction='in', width=2.5, length=12, zorder=1, right=True)
    this_ax.xaxis.set_tick_params(which='minor', direction='in', width=1.5, length=6, zorder=1, top=True)
    this_ax.yaxis.set_tick_params(which='minor', direction='in', width=1.5, length=6, zorder=1, right=True)
    dx = -3/72
    dy = -3/72
    y_offset = matplotlib.transforms.ScaledTranslation(0, dy, this_fig.dpi_scale_trans)
    x_offset = matplotlib.transforms.ScaledTranslation(dx, 0, this_fig.dpi_scale_trans)

    for label in this_ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + y_offset)

    for label in this_ax.yaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + x_offset)

        