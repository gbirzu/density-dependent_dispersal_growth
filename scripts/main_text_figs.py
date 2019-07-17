import os
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import analytic_tools as analytic
import analysis_tools as analysis
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


phase_plot_file = '../data/data_2d_linear-migration_cooperative-growth.txt'
migration_file = '../data/data_r0-0.001_m0-0.05_logistic.npy'
input_path = '../data/'
output_path = '../figures/'


def make_fig1_panels( ne = None, save=True ):
    # ---------------------------------------------------------------------------------- #
    profile_fname = '../data/fluctuations/N10000/N10000_r0.001_m0.05_A0.0_B0.0/profile/profile_N10000_r0.001_m0.05_A0.0_B0.0_run85.txt'
    t_array, n_array, f_array, f_mean, f_std = analysis.read_fraction_distribution( profile_fname, n_cutoff=0 )
    x_array = np.arange( n_array.shape[ 1 ] )

    t_index_list = [ 0, 2, 20, 60 ]
    delta_index = 25
    n_max = np.max( n_array )

    for i in range( 4 ):
        fig = plt.figure( figsize=analysis.cm2inch( 8.7, 6.0 ) )
        ax = fig.add_subplot( 111 )

        t_index = t_index_list[ i ]
        rhot_array = n_array[ t_index ] / n_max
        ft_array = f_array[ t_index ]
        rho1t_array = ft_array * rhot_array
        i_index = 47
        f_index = 202
        x_plot = x_array[ i_index:f_index ]
        t = t_array[ t_index ]
        if ne is not None:
            t /= ne

        print(t)

        # Setup axes
        ax.set_xticks( [ 50, 100, 150, 200 ] )
        ax.set_yticks( [ 0, 0.5, 1 ] )
        ax.set_xticklabels( [] )
        ax.set_yticklabels( [] )

        ax.set_xlim( [ x_plot[ 0 ], x_plot[ -1 ] ] )
        ax.set_ylim( [ 0.0, 1.2 ] )

        # Plot data
        ax.plot( x_plot, rhot_array[ i_index:f_index ], lw=2, c='k' )
        ax.plot( x_plot, rho1t_array[ i_index:f_index ], lw=1.0, ls='-', c='k' )
        ax.fill_between( x_plot, rho1t_array[ i_index:f_index ], rhot_array[ i_index:f_index ], facecolor='r', alpha=0.5 )
        ax.fill_between( x_plot, np.zeros( f_index - i_index ), rho1t_array[ i_index:f_index ], facecolor='g', alpha=0.5 )

        if save == True:
            plt.savefig( output_path + 'fig1/fig1_panel_pulled' + str( i + 1 ) + '.pdf' )
    # ---------------------------------------------------------------------------------- #

    fig = plt.figure( figsize=analysis.cm2inch( 8.7, 6.0 ) )
    ax = fig.add_subplot( 111 )

    t_index = 20
    rhot_array = n_array[ t_index ] / n_max
    ft_array = f_array[ 2 ]
    rho1t_array = ft_array * rhot_array
    i_index = 47
    f_index = 202
    x_plot = x_array[ i_index:f_index ]

    t = t_array[ t_index ]
    if ne is not None:
        t /= ne

    print(t)

    # Setup axes
    ax.set_xticks( [ 50, 100, 150, 200 ] )
    ax.set_yticks( [ 0, 0.5, 1 ] )
    ax.set_xticklabels( [] )
    ax.set_yticklabels( [] )

    ax.set_xlim( [ x_plot[ 0 ], x_plot[ -1 ] ] )
    ax.set_ylim( [ 0.0, 1.2 ] )

    # Plot data
    ax.plot( x_plot, rhot_array[ i_index:f_index ], lw=2, c='k' )
    ax.plot( x_plot, rho1t_array[ i_index:f_index ], lw=1.0, ls='-', c='k' )
    ax.fill_between( x_plot, rho1t_array[ i_index:f_index ], rhot_array[ i_index:f_index ], facecolor='r', alpha=0.5 )
    ax.fill_between( x_plot, np.zeros( f_index - i_index ), rho1t_array[ i_index:f_index ], facecolor='g', alpha=0.5 )

    if save == True:
        plt.savefig( output_path + 'fig1/fig1_panel_fluctuations' + str( i + 1 ) + '.pdf' )


    fig = plt.figure( figsize=analysis.cm2inch( 8.7, 6.0 ) )
    ax = fig.add_subplot( 111 )

    t_index = 2
    rhot_array = analytic.theory_profile(x_array, 165, 0.15)
    ft_array = f_array[ t_index ]
    rho1t_array = ft_array * rhot_array
    i_index = 47
    f_index = 202
    x_plot = x_array[ i_index:f_index ]

    t = t_array[ t_index ]
    if ne is not None:
        t /= ne

    print(t)

    # Setup axes
    ax.set_xticks( [ 50, 100, 150, 200 ] )
    ax.set_yticks( [ 0, 0.5, 1 ] )
    ax.set_xticklabels( [] )
    ax.set_yticklabels( [] )

    ax.set_xlim( [ x_plot[ 0 ], x_plot[ -1 ] ] )
    ax.set_ylim( [ 0.0, 1.2 ] )

    # Plot data
    ax.plot( x_plot, rhot_array[ i_index:f_index ], lw=2, c='k' )
    ax.plot( x_plot, rho1t_array[ i_index:f_index ], lw=1.0, ls='-', c='k' )
    ax.fill_between( x_plot, rho1t_array[ i_index:f_index ], rhot_array[ i_index:f_index ], facecolor='r', alpha=0.5 )
    ax.fill_between( x_plot, np.zeros( f_index - i_index ), rho1t_array[ i_index:f_index ], facecolor='g', alpha=0.5 )

    if save == True:
        plt.savefig( output_path + 'fig1/fig1_panel_theory' + str( i + 1 ) + '.pdf' )


def cooperation( save=True ):
    matplotlib.rcParams.update( { 'font.size':9 } )
    label_size = 12
    fig = plt.figure( figsize=analysis.cm2inch( 17.8, 8.5 ) )

    # Plotting tools
    symbol = [ 'o', 's', 'D' ]
    colors = [ '0.6', '0.3', '0.0' ]
    n_array = [ 10000, 100000, 1000000 ]
    n_dict = { 10000:'$10^4$', 100000:'$10^5$', 1000000:'$10^6$' }
    # ---------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------- #
    ax = fig.add_subplot( 121 )
    ax.set_xlabel( 'time, t', fontsize=label_size, fontweight='bold', color='k' )
    ax.set_ylabel( 'heterozygosity, H', fontsize=label_size, fontweight='bold', color='k' )
    ax.set_yscale( 'log' )
    ax.ticklabel_format( style='sci', scilimit=(-2,2), axis='x' )
    ax.set_xticks( [ 0, 20000, 40000, 60000, 80000 ] )

    # Load data
    fname_list = [ 'hetero_N10000_r0.001_m0.05_A0.0_B0.0_avg.npy' ]

    # Set plot parameters
    symbol_size = 40
    number_of_points = 30

    for fname in fname_list:
        het_data = np.load( input_path + 'cooperation/r0.001_m0.05/' + fname )
        time_array = het_data[ 0 ]
        het_array = het_data[ 1 ]
        survival_array = het_data[ 2 ]

        # Keep only nonzero data
        time_array = time_array[ np.where( het_array > 1E-3 )[ 0 ] ]
        survival_array = survival_array[ np.where( het_array > 1E-3 )[ 0 ] ]
        het_array = het_array[ np.where( het_array > 1E-3 )[ 0 ] ]

        # Filter intermediate points
        reduced_het_indices = ( len( het_array )/number_of_points )*np.arange( number_of_points )
        het_plot = np.array( [ het_array[ i ] for i in reduced_het_indices ] )
        x_plot = np.array( [ time_array[ i ] for i in reduced_het_indices ] )

        # Fit H(t) to exponential
        xh, het_fit, coeffs, r_sq = analysis.fit_heterozygosity( time_array, het_array, survival_array )
        fit = np.poly1d( coeffs )
        t_begin = time_array[ np.where( het_array < 0.4 )[ 0 ][ 0 ] ]
        t_end = 1.1*time_array[ reduced_het_indices[ -1 ] ]
        x_fit = [ t_begin, t_end ] #choose range for plotting fit
        est = np.exp( fit( x_fit ) )

        # Plot data
        ax.scatter( x_plot, het_plot, s=symbol_size, edgecolor='none', facecolor='gray', label='simulations' )
        ax.plot( x_fit, est, c='k', lw=1, label='fit' )
    ax.legend( loc='upper right' )
    # ---------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------- #
    # Get data
    raw_data = np.load( '../data/cooperation/r0.001_m0.05/neff_data_phaseplot.npy' )
    allee_data = []
    diffusion_data = []
    for elem in raw_data:
        for row in elem:
            if row[ 2 ] == 0.0:
                allee_data.append( row )
            elif row[ 3 ] == 0.0:
                diffusion_data.append( row )
    allee_data = np.array( allee_data )
    diffusion_data = np.array( diffusion_data )

    # Define plot parameters
    A_min = -0.6
    A_max = 2.9
    Ne_min = -5E4
    Ne_max = 3.2E6

    # Setup axis
    ax = fig.add_subplot( 122 )
    ax.set_xlim( A_min, A_max )
    ax.set_ylim( Ne_min, Ne_max )
    ax.set_xticks( [ -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 ] )
    ax.set_xticklabels( [ '-0.5', '0', '0.5', '1', '1.5', '2', '2.5' ] )
    ax.set_yticks( [ 0, 1E6, 2E6, 3E6 ] )
    ax.set_yticklabels( [ '$\mathbf{0}$', '$\mathbf{1 \cdot 10^6}$', '$\mathbf{2 \cdot 10^6}$', '$\mathbf{3 \cdot 10^6}$' ] )
    ax.set_xlabel( "dispersal cooperativity, $\mathbf{A_1}$", fontweight="bold", fontsize=label_size )
    ax.set_ylabel( "effective population size, $\mathbf{N_e}$", fontweight="bold", fontsize=label_size )

    # Plot data
    diffusion_data = diffusion_data[ diffusion_data[ :, 2 ].argsort() ] # Sort by value of B
    for i, n in enumerate( n_array[ 2:3 ] ):
        indices = np.where( diffusion_data.T[ 0 ] == n )[ 0 ] # Find entries with N population
        a_array = diffusion_data.T[ 2 ][ indices ]
        ne_array = diffusion_data.T[ 1 ][ indices ]
        f = interpolate.CubicSpline(a_array, ne_array)
        x_new = np.linspace(a_array[0], a_array[-1], num=50, endpoint=True)
        ax.plot( x_new, f(x_new), ls='-', lw=1, c='k' )
        ax.scatter( a_array, ne_array, marker=symbol[ i ], s=25.0, c='k', facecolor='k', edgecolor='none', label='simulations' )
    # ---------------------------------------------------------------------------------- #

    ax.legend( loc='upper left' )

    plt.tight_layout()
    if save==True:
        plt.savefig( output_path + "Fig2_cooperation.pdf" )


def make_cooperation_panels( save=True ):
    matplotlib.rcParams.update( { 'font.size':8 } )
    label_size = 12
    fig = plt.figure( figsize=analysis.cm2inch( 6.8, 5.4 ) )

    # Plotting tools
    symbol = [ 'o', 's', 'D' ]
    colors = [ '0.6', '0.3', '0.0' ]
    n_array = [ 10000, 100000, 1000000 ]
    n_dict = { 10000:'$10^4$', 100000:'$10^5$', 1000000:'$10^6$' }
    # ---------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------- #
    ax = fig.add_subplot( 111 )
    ax.set_yscale( 'log' )
    ax.ticklabel_format( style='sci', scilimit=(-2,2), axis='x' )
    ax.set_xticks( [ 0, 2, 4, 6 ] )
    ax.set_xticklabels([])
    ax.set_yticks([1E-3, 1E-2])
    ax.set_yticklabels([])

    # Load data
    fname_list = [ 'hetero_N10000_r0.001_m0.05_A0.0_B0.0_avg.npy' ]

    # Set plot parameters
    symbol_size = 40
    number_of_points = 30
    ne=0

    for fname in fname_list:
        het_data = np.load( input_path + 'cooperation/r0.001_m0.05/' + fname )
        time_array = het_data[ 0 ]
        het_array = het_data[ 1 ]
        survival_array = het_data[ 2 ]

        # Keep only nonzero data
        time_array = time_array[ np.where( het_array > 1E-3 )[ 0 ] ]
        survival_array = survival_array[ np.where( het_array > 1E-3 )[ 0 ] ]
        het_array = het_array[ np.where( het_array > 1E-3 )[ 0 ] ]

        # Filter intermediate points
        reduced_het_indices = ( len( het_array )/number_of_points )*np.arange( number_of_points )
        het_plot = np.array( [ het_array[ i ] for i in reduced_het_indices ] )
        x_plot = np.array( [ time_array[ i ] for i in reduced_het_indices ] )

        # Fit H(t) to exponential
        xh, het_fit, coeffs, r_sq = analysis.fit_heterozygosity( time_array, het_array, survival_array )
        ne = -1 / coeffs[ 0 ]
        time_array /= ne
        x_plot /= ne
        fit = np.poly1d( [-1, coeffs[1]] )
        t_begin = time_array[ np.where( het_array < 0.4 )[ 0 ][ 0 ] ]
        t_end = 1.1*time_array[ reduced_het_indices[ -1 ] ]
        x_fit = np.array( [ t_begin, t_end ] ) #choose range for plotting fit
        est = np.exp( fit( x_fit ) )

        # Plot data
        ax.scatter( x_plot, het_plot, s=symbol_size, edgecolor='none', facecolor='gray', label='simulations' )
        ax.plot( x_fit, est, c='k', lw=1, label='fit' )
    # ---------------------------------------------------------------------------------- #

    plt.tight_layout()
    if save==True:
        plt.savefig( output_path + "fig1/cooperation_panel-a_timescale-ne.pdf" )

    return ne


def theory_summary( save=True ):
    font = {'family' : 'sans-serif', 'serif' : 'Helvetica Neue', 'weight' : 'bold', 'size' : 10}
    matplotlib.rc( 'font', **font )

    # Set plot parameters
    label_size = 11

    rho = np.arange( 0, 1.01, 0.01 )

    fig = plt.figure( figsize=analysis.cm2inch( 18.4, 9.0 ) )

    # ---------------------------------------------------------------------------------- #
    BMin = 0
    BMax = 6.0
    B_array = np.arange( BMin, BMax, 0.001 )
    B_pulled_array = np.arange( BMin, 2.0, 0.001 )
    B_semipushed_array = np.arange( 2.0, 4.0, 0.001 )
    B_fullypushed_array = np.arange( 4.0, BMax, 0.001 )
    v_array = np.array( [ analytic.vA( 1.0, 1.0, B ) for B in B_array ] )
    v_fisher = min( v_array )

    v_min = 0.95*min( v_array )
    v_max = 1.05*max( v_array )

    y_pp_transition = np.arange( v_min, v_max, 0.001 )
    x_pp_transition = -0.5*np.ones( len( y_pp_transition ) )

    ax = fig.add_subplot( 121 )
    ax.set_xlim( [ BMin, BMax ] )
    ax.set_ylim( [ v_min, v_max ] )
    ax.set_yticks( [] )
    ax.set_xlabel( 'cooperativity, A or B', fontsize=label_size, fontweight='bold' )
    ax.set_ylabel( 'normalized velocity, $\mathbf{\\nu}$', fontsize=label_size, fontweight='bold' )

    ax.set_xticks( [] )
    ax.tick_params( axis='both', which='major', labelsize=9 )
    ax.tick_params( axis='both', which='minor', labelsize=9 )

    ax.fill_between( B_pulled_array, v_min, v_max, facecolor='lightsalmon', alpha=0.7 )
    ax.fill_between( B_semipushed_array, v_min, v_max, facecolor='lightgreen', alpha=0.7 )
    ax.fill_between( B_fullypushed_array, v_min, v_max, facecolor='lightskyblue', alpha=0.7 )

    ax.plot( B_array, v_array, ls='-', lw=2, c='k' )
    ax.plot( B_array, 2.*np.ones( len( B_array ) ), ls='--', lw=2, c='k' )
    ax.plot( B_array, ( 3/( 2*np.sqrt(2) ) )*v_fisher*np.ones( len( B_array ) ), ls='--', lw=2, c='k' )
    ax.text( 1.00, 1.015*v_max, 'pulled', horizontalalignment="center", fontsize=11, fontweight='bold', color='lightsalmon' )
    ax.text( 3.00, 1.010*v_max, 'semi-' + '\npushed', horizontalalignment="center", fontsize=11, fontweight='bold', color='lightgreen' )
    ax.text( 5.00, 1.010*v_max, 'fully-' + '\npushed', horizontalalignment="center", fontsize=11, fontweight='bold', color='lightskyblue' )
    # ---------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------- #
    ax = fig.add_subplot( 122 )

    # Define parameters
    ax.set_xlim( [ 2E2, 6E9 ] )
    ax.set_ylim( [ 2E2, 1E9 ] )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xticks( [ 1E3, 1E5, 1E7, 1E9 ] )
    ax.set_yticks( [ 1E3, 1E5, 1E7, 1E9 ] )
    ax.set_xlabel( 'population size, $\mathbf{N}$', fontsize=label_size, fontweight='bold' )
    ax.set_ylabel( 'effective population size, $\mathbf{N_e}$', fontsize=label_size, fontweight='bold' )

    #ax.set_xticks( [] )
    ax.tick_params( axis='both', which='major', labelsize=9 )
    ax.tick_params( axis='both', which='minor', labelsize=9 )

    # Define data
    n_max = 1E9
    n_array = np.logspace( np.log10( 5E2 ), np.log10( n_max/20 ), base=10 )
    ne_pulled = np.log( n_array )**3
    ne_data_min = 5E2
    ne_pulled = ne_data_min*ne_pulled/min( ne_pulled )
    ne_fullypushed = n_array
    ne_semipushed = np.power( n_array, 0.6 )
    ne_semipushed = ( ( ne_fullypushed[ 0 ] + ne_pulled[ 0 ] )/2 )*ne_semipushed/min( ne_semipushed )

    exponent_array = [ 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
    ne_range = ne_fullypushed[ 0 ] - ne_pulled[ 0 ]
    delta_ne = ne_range/( len( exponent_array ) + 1 )

    ax.plot( n_array, ne_pulled, ls='-', lw=2, c='r' )
    ax.plot( n_array, ne_fullypushed, ls='-', lw=2, c='b' )

    for i, exponent in enumerate( exponent_array ):
        ne_semipushed = np.power( n_array, exponent )

        # Rescale Ne for clarity
        ne_semipushed_min = ne_data_min + ( i + 1 )*delta_ne
        ne_semipushed = ne_semipushed_min*ne_semipushed/min( ne_semipushed )
        ax.plot( n_array, ne_semipushed, ls='--', lw=2, c='g' )

    ax.text( 1E4, 6.0E6, 'fully-pushed', color='b', rotation=42, fontsize=11 )
    ax.text( 7E7, 5E5, 'semi-\npushed', color='g', rotation=0, fontsize=11 )
    ax.text( 1E5, 1.6E3, 'pulled', color='r', rotation=8, fontsize=11 )
    # ---------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------- #
    plt.tight_layout( pad=3.0, w_pad=1.0, h_pad=2.5 )
    if save == True:
        plt.savefig( output_path + 'Fig3_theory_summary.pdf' )


def diversity( theory=True, save=True ):
    font = {'family' : 'sans-serif', 'serif' : 'Helvetica Neue', 'weight' : 'bold', 'size' : 10}
    matplotlib.rc( 'font', **font )

    axis_fontsize = 12

    het_data = np.load( input_path + 'hetero_comparison_data.npy' )

    fig = plt.figure( figsize=analysis.cm2inch( 17.8, 9.0 ) )

    ax1 = fig.add_subplot( 121 )
    ax1.set_xlim( -0.03, 1.05 )
    ax1.set_ylim( 0.0, 0.21 )
    ax1.set_xticks( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
    ax1.set_yticks( [ 0.0, 0.05, 0.1, 0.15, 0.2 ] )
    ax1.set_xlabel( 'population density, $\mathbf{n/N}$', fontweight='bold', fontsize=axis_fontsize )
    ax1.set_ylabel( 'dispersal rate, $\mathbf{D( n )}$', fontweight='bold', fontsize=axis_fontsize )
    m0 = 0.1
    A = 1
    rho = np.arange( 0, 1.1, 0.10 )
    rho_line = np.linspace( 0, 1.0, 1000 )
    ax1.scatter( rho, m0*( 1. + A*rho ), lw=1, marker='o', facecolor='none', edgecolor='k' )
    ax1.plot( rho_line, m0*( 1. + A*rho_line ), lw=0.5, c='k' )
    ax1.scatter( rho, m0*( 1. + A*np.sqrt( rho ) ), lw=1, marker='^', facecolor='none', edgecolor='k' )
    ax1.plot( rho_line, m0*( 1. + A*np.sqrt( rho_line ) ), lw=0.5, c='k' )
    ax1.scatter( rho, m0*( 1. + A*rho**2 ), lw=1, marker='s', facecolor='none', edgecolor='k' )
    ax1.plot( rho_line, m0*( 1. + A*rho_line**2 ), lw=0.5, c='k' )
    ax1.plot( rho, m0*np.ones( len( rho ) ), ls="--", lw=1, c='k' )
    ax1.plot( rho, 2*m0*np.ones( len( rho ) ), ls="--", lw=1, c='k' )
    # ---------------------------------------------------------------------------------- #

    # ---------------------------------------------------------------------------------- #
    ax = fig.add_subplot( 122 )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xlabel( 'population size, $\mathbf{N}$', fontweight='bold', fontsize=axis_fontsize )
    ax.set_ylabel( 'effective population size, $\mathbf{N_e}$', fontweight='bold', fontsize=axis_fontsize )
    ax.set_xticks( [ 1E3, 1E4, 1E5, 1E6, 1E7 ] )
    ax.set_yticks( [ 1E3, 1E4, 1E5, 1E6, 1E7, 1E8 ] )
    ax.set_xlim( [ 2.0E3, 2.5E7 ] )
    ax.set_ylim( [ 1.5E3, 2.5E8 ] )

    ax_inset = inset_axes( ax, width=1.0, height=0.9, loc=2, borderpad=0.25 )
    ax_inset.set_xscale( 'log' )
    ax_inset.set_yscale( 'log' )
    ax_inset.set_xlabel( '$\mathbf{N}$', fontweight='bold', fontsize=8, labelpad=0 )
    ax_inset.yaxis.set_label_position( 'right' )
    ax_inset.set_ylabel( '$\mathbf{N_e}$', fontweight='bold', fontsize=8, rotation=0, labelpad=6 )
    ax_inset.set_xticks( [ 1E4, 1E5, 1E6 ] )
    ax_inset.set_yticks( [ 1E4, 1E5, 1E6, 1E7 ] )
    ax_inset.tick_params( axis='y', which='both', left='off', labelleft='off', right='on', labelright='on', labelsize=7 )
    ax_inset.tick_params( axis='x', labelsize=7 )
    ax_inset.set_xlim( [ 3.0E3, 8.0E6 ] )
    ax_inset.set_ylim( [ 2.0E3, 2.8E7 ] )

    markers = [ '^', 'o', 's' ]
    labels = [ '$m_0( 1 + A_{1/2}\sqrt{\\rho} )$', '$m_0( 1 + A_1\\rho )$', '$m_0( 1 + A_2\\rho^2 )$' ]
    colors = [ 'red', 'green', 'blue' ]
    model_names = [ 'sqrt', 'linear', 'quadratic' ]

    sqrt_data = np.load( input_path + 'scaling_data/sqrt_model/neff_data_phaseplot.npy' )
    linear_data = np.load( input_path + 'scaling_data/linear_model/neff_data_phaseplot.npy' )
    quadratic_data = np.load( input_path + 'scaling_data/quadratic_model/neff_data_phaseplot.npy' )
    full_data = [ sqrt_data, linear_data, quadratic_data ]
    df_params = pd.read_csv( '../data/fig4_parameters.csv', sep=',' )

    nu_list = [ 1.0, 1.03, 1.5 ]
    for i, nu in enumerate( nu_list ):
        for j, model in enumerate( full_data ):
            # Get data for nu
            dispersal_coop = df_params.loc[ df_params[ 'nu' ] == nu, 'm1/m0_' + model_names[ j ] ].values[ 0 ]
            neff_data = full_data[ j ][ np.where( full_data[ j ][ :, 0, 2 ] == dispersal_coop )[ 0 ][ 0 ] ]
            neff_data = neff_data[ np.argsort( neff_data[ :, 0 ] ) ]
            neff_data = neff_data[ np.where( neff_data[ :, 0 ] >= 10000 )[ 0 ] ]

            # Get N_array and Ne_array
            neff_array = neff_data[ :, 1 ]
            n_array = neff_data[ :, 0 ]

            i_min = 2
            i_fit = 3
            # Normalize to sqrt model
            neff_collapse_array = neff_array
            if i == 0:
                neff_collapse_array *= 0.5

            if j == 0:
                neff_scaling = neff_collapse_array
            else:
                neff_collapse_array = neff_collapse_array * neff_scaling[ i_fit ] / neff_collapse_array[ i_fit ]


            # Fit slope and plot
            slope, intercept, r_value, p_value, std_err = stats.linregress( np.log( n_array[ i_min: ] ), np.log( neff_array[ i_min: ] ) )
            n_fit = np.array( [ 0.5*min( n_array ), n_array[ i_fit ], 2.5*max( n_array ) ] )

            if theory == True:
                theory_exponent = analytic.exact_exponent( nu )
                log_weight = int( theory_exponent == 0 )
                neff_fit = ( 1 - log_weight ) * n_fit ** theory_exponent + log_weight * ( np.log( n_fit ) ) ** 3
                neff_fit_collapse = neff_fit * neff_collapse_array[ i_fit ] / neff_fit[ 1 ]
                neff_fit = neff_fit * neff_array[ i_fit ] / neff_fit[ 1 ]
            else:
                neff_fit = np.exp( slope*np.log( n_fit ) + intercept )

            if i == 0:
                ax.plot( n_fit, neff_fit, lw=1, c=colors[ i ], zorder=1 )
                ax.scatter( n_array, neff_array, marker=markers[ j ], edgecolor=colors[ i ], facecolor='none', s=30, label=labels[ j ], lw=1, zorder=2 )
                ax_inset.plot( n_fit, neff_fit_collapse, lw=1, c=colors[ i ], zorder=1 )
                ax_inset.scatter( n_array, neff_collapse_array, marker=markers[ j ], edgecolor=colors[ i ], facecolor='none', s=30, label=labels[ j ], lw=1, zorder=2 )
            else:
                ax.plot( n_fit, neff_fit, lw=1, c=colors[ i ], zorder=1 )
                ax.scatter( n_array, neff_array, marker=markers[ j ], edgecolor=colors[ i ], facecolor='none', s=30, lw=1, zorder=2 )
                ax_inset.plot( n_fit, neff_fit_collapse, lw=1, c=colors[ i ], zorder=1 )
                ax_inset.scatter( n_array, neff_collapse_array, marker=markers[ j ], edgecolor=colors[ i ], facecolor='none', s=30, label=labels[ j ], lw=1, zorder=2 )


    ax.text( 3E6, 1.8E7, 'fully-\npushed', color='b', fontsize=10 )
    ax.text( 3E6, 9E5, 'semi-\npushed', color='g', fontsize=10 )
    ax.text( 3E6, 5.0E4, 'pulled', color='r', fontsize=10 )

    # Create legend
    legend_elements = [ Line2D( [ 0 ], [ 0 ], marker='^', color='w', label='$m_0 \left( 1 + A_{1/2}\sqrt{n/N} \,\, \\right)$', markeredgecolor='k', markerfacecolor='none', markersize=6, markeredgewidth=1 ),
                            Line2D( [ 0 ], [ 0 ], marker='o', color='w', label='$m_0 \left( 1 + A_1 n/N \\right)$', markeredgecolor='k', markerfacecolor='none', markersize=6, markeredgewidth=1 ),
                            Line2D( [ 0 ], [ 0 ], marker='s', color='w', label='$m_0 \left[ 1 + A_2\left( n/N \\right)^2 \\right]$', markeredgecolor='k', markerfacecolor='none', markersize=6, markeredgewidth=1 ) ]

    legend_properties={'weight':'bold', 'size':9}
    ax1.legend( handles=legend_elements, loc='lower right', prop=legend_properties, scatterpoints=1, handletextpad=0.2 )

    plt.tight_layout( pad=0.5 )
    plt.savefig( output_path + 'Fig4_diversity.pdf' )


def universality( export_parameters=False, save=True ):
    font = { 'family':'sans-serif', 'serif':'Helvetica Neue', 'weight':'bold', 'size':9 }
    matplotlib.rc( 'font', **font )

    # Set plot parameters
    label_size = 10

    # Define constants
    r0 = 0.001
    m0 = 0.05
    v_F = 2.*np.sqrt( r0*m0/2. )
    delta_s = 3./( 2*np.sqrt( 2 ) ) - 1.
    A_min = 0.0
    A_max = 2.0

    # Load data
    data_array = np.loadtxt( phase_plot_file, delimiter=',' )

    x = np.sort( np.unique( data_array.T[ 0 ] ) )
    y = np.sort( np.unique( data_array.T[ 1 ] ) )
    v_array = data_array.T[ 2 ]
    [ x_matrix, y_matrix, v_matrix ], v_values, v_plot = analysis.plot_matrix( x, y, v_array, v_F )

    v_fisher = min( v_array )
    v_semi = ( 1 + delta_s/2.0 )*v_fisher
    v_fully = 1.30*v_fisher # Minimum velocity to avoid crossovers
    v_contour = [ v_semi, 1.15*v_fisher, v_fully, 1.45*v_fisher, 1.6*v_fisher ]
    epsilon = 1E-4

    # Choose parameters for simulations
    all_params = [ np.array( [ 0.0, 0.5, 1 ] ), np.array( [ 0.125, 0.25, 1 ] ), np.array( [ 0.25, 0.0, 1 ] ) ]

    f_phaseplot = open( 'ab_parameters.txt', 'w' )
    # Semi-pushed waves
    params_semi = data_array[ np.where( abs( v_array - v_semi )/v_semi < epsilon )[ 0 ] ]

    for elem in params_semi:
        f_phaseplot.write( str( elem[ 0 ] ) + ',' + str( elem[ 1 ] ) + '\n' )

    # Find linear approximation
    sorted_params = params_semi[ np.argsort( params_semi[ :, 0 ] ) ]
    [ x0, y0 ] = sorted_params[ 0 ][ :2 ]
    [ x1, y1 ] = sorted_params[ -1 ][ :2 ]
    xm = ( x0 + x1 )/2.
    ym = ( y0 + y1 )/2.
    i_x = np.where( abs( x - xm ) < 0.003 )[ 0 ][ 0 ]
    i_y = np.where( abs( y - ym ) < 0.003 )[ 0 ][ 0 ]
    v_linear = v_matrix[ i_x, i_y ]

    # Pick plot parameters
    params_semi = analysis.pick( params_semi, 5, int( len( params_semi )*0.1 ) )
    for params in params_semi:
        all_params.append( params )

    # Fully-pushed waves
    params_fully = data_array[ np.where( abs( v_array - v_fully )/v_fully < epsilon )[ 0 ] ]
    for elem in params_fully:
        f_phaseplot.write( str( elem[ 0 ] ) + ',' + str( elem[ 1 ] ) + '\n' )
    f_phaseplot.close()
    params_fully = analysis.pick( params_fully, 7, int( len( params_fully )*0.1 ) )
    for params in params_fully:
        all_params.append( params )
    all_params = np.array( all_params )

    # Export simulations parameters
    if export_parameters == True:
        f_parameters = open( "2d_parameters.csv", "w" )
        for params in all_params:
            f_parameters.write( str( params[ 0 ] ) + "," + str( params[ 1 ] ) + "\n" )
        f_parameters.close()

    fig = plt.figure( figsize=analysis.cm2inch( 17.8, 9.0 ) )
    cmap = {v_values[ 0 ]:analysis.colors_dict[ 'pulled' ], v_values[ 1 ]:analysis.colors_dict[ 'semi-pushed' ], v_values[ 2 ]:analysis.colors_dict[ 'fully-pushed' ]}
    cm = matplotlib.colors.LinearSegmentedColormap.from_list( 'phase plot colors', analysis.colors, N=3 )

    ax = fig.add_subplot( 121 )
    ax.set_xlim( 0, 2.5 )
    ax.set_ylim( 0, 5 )
    ax.set_xlabel( "dispersal cooperativity, $\mathbf{A_1}$", fontweight="bold", fontsize=label_size )
    ax.set_ylabel( "growth cooperativity, $\mathbf{B}$", fontweight="bold", fontsize=label_size )

    ax.imshow( v_plot, cmap=cm, origin='lower', extent=[ x[ 0 ], x[ -1 ], y[ 0 ], y[ -1 ] ], aspect=0.5 )
    ax.contour( x_matrix.T, y_matrix.T, v_matrix, levels=v_contour, colors='k', linewidths=0.5 )
    ax.scatter( all_params[ :, 0 ], all_params[ :, 1 ], s=20, lw=2, facecolor='k', marker='x' )

    ax = fig.add_subplot( 122 )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xlim( [ 3E3, 3E6 ] )
    ax.set_xlim( [ 3E3, 3E6 ] )
    ax.set_xlabel( "population size, N", fontweight="bold", fontsize=label_size )
    ax.set_ylabel( "effective population size, $\mathbf{N_e}$", fontweight="bold", fontsize=label_size )

    neff_data = np.load( '../data/neff_data_phaseplot.npy' )
    for experiment in neff_data:
        A = experiment[ 0 ][ 2 ]
        B = experiment[ 0 ][ 3 ]

        n_array = experiment.T[ 0 ]
        neff_array = experiment.T[ 1 ]

        if 2*A + B <= 2.0:
            #clr = 'r'
            clr = 'lightsalmon'
            marker_clr = 'r'
        elif 2*A + B < 4.0:
            #clr = 'g'
            clr = 'lightgreen'
            marker_clr = 'g'
        else:
            #clr = 'b'
            clr = 'lightskyblue'
            marker_clr = 'b'

        n_fit = np.array( [ 0.5*min( n_array ), n_array[ len( n_array )/2 ], 1.5*max( n_array ) ] )
        slope, intercept, r_value, p_value, std_err = stats.linregress( np.log( n_array ), np.log( neff_array ) )
        ax.plot( n_fit, np.exp( intercept + slope*np.log( n_fit ) ), lw=1, c=clr, zorder=1 )
        ax.scatter( n_array, neff_array, marker='x', facecolor=marker_clr, edgecolor='none', s=20, zorder=2 )

        print A, B, slope

    alpha_theory = analytic.exact_exponent( v_semi/v_fisher )
    alpha_linear = analytic.exact_exponent( v_linear/v_fisher )
    N_test = 1E9
    print alpha_theory, N_test**alpha_theory
    print alpha_linear, N_test**alpha_linear

    plt.tight_layout( pad=3.0, w_pad=1.0 )

    # Create legend
    legend_elements = [ Line2D( [ 0 ], [ 0 ], marker='o', color='w', label='pulled', markerfacecolor='lightsalmon', markersize=10 ),
                        Line2D( [ 0 ], [ 0 ], marker='o', color='w', label='semi-pushed', markerfacecolor='lightgreen', markersize=10 ),
                        Line2D( [ 0 ], [ 0 ], marker='o', color='w', label='fully-pushed', markerfacecolor='lightskyblue', markersize=10 ) ]
    ax.legend( handles=legend_elements, loc='center', bbox_to_anchor=( -0.20, 1.1 ), ncol=3, borderpad=0.5, edgecolor='k', shadow=False )

    if save == True:
        plt.savefig( output_path + 'Fig5_universality.pdf', dpi=500 )


if __name__ == '__main__':
    font = {'family' : 'sans-serif', 'serif' : 'Helvetica Neue', 'weight' : 'bold', 'size' : 8}
    matplotlib.rc( 'font', **font )

    ne = make_cooperation_panels( save=True )
    make_fig1_panels( ne = ne, save=True )
    cooperation( save=True )
    theory_summary( save=True )
    diversity( save=True )
    universality( export_parameters=False, save=True )

    plt.show()
