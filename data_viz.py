import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def boxplot(data, meta, x_label, y_label, title, out_file):
    """plot boxplot for input parallel array and save the result as a png file
    """
    # settings for plotting
    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    fig, ax1 = plt.subplots(figsize=(8, 6))

    ax1.boxplot(data)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_title(title)
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    # set tick labels and export the result
    ax1.set_xticklabels(meta, rotation=90, fontsize=8)
    fig.savefig(out_file, bbox_inches='tight')
