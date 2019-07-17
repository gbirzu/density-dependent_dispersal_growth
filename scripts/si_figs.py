import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
import scipy.stats as stats

data_path = '../data/het_average.dat'
output_dir = '../figures/'

# Configure matplotlib environment
helvetica_scale_factor = 0.92 # rescale Helvetica to other fonts of same size
mpl.rcParams['font.size'] = 10 * helvetica_scale_factor
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Helvetica Neue'
mpl.rcParams['axes.titlesize'] = 12 * helvetica_scale_factor

single_col_width = 3.43 # = 8.7 cm
double_col_width = 7.01 # = 17.8 cm


def plot_het_comparison(het_averages):
    time = het_averages['time']
    het_global = het_averages['global']
    het_local = het_averages['local']

    fig = plt.figure(figsize=(single_col_width, single_col_width))
    ax = fig.add_subplot(111)
    ax.set_xlabel('time, t', fontweight='bold')
    ax.set_ylabel('heterozygosity, H', fontweight='bold')
    ax.set_yscale('log')
    ax.plot(time, het_global, ls='-', lw=2, c='k')
    ax.plot(time, het_local, ls='', marker='o', markevery=5, markersize=5, markeredgecolor='r', markerfacecolor='none')

    plt.tight_layout()
    plt.savefig(output_dir + 'het_comparison.pdf')


def fit_Ne(het_averages, averaging='global'):
    time = het_averages['time']
    het = het_averages[averaging]

    slope, intercept, rvalue, pvalue, stderr = stats.linregress(time, np.log(het))
    return 1 / abs(slope)


if __name__ == '__main__':
    with open(data_path, 'rb') as f_in:
        het_averages = pickle.load(f_in)
    plot_het_comparison(het_averages)
    ne_global = fit_Ne(het_averages, averaging='global')
    ne_local = fit_Ne(het_averages, averaging='local')
    print('Ne (global averaging): ', ne_global)
    print('Ne (local averaging): ', ne_local)
    print('Ne difference: ', 100 * (ne_global - ne_local) / ne_global, '%')
