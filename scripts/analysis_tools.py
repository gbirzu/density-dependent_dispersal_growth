import glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import analytic_tools as analytic
from scipy import stats
from scipy import optimize

input_path = '../data/'
output_plots = '../plots/'
output_files = '../data/'

# Set up colors dict
#colors_dict = {'lightsalmon':( 1., 0.627, 0.478 ), 'lightskyblue':( 0.529, 0.808, 0.98 ), 'lightgreen':( 0.565, 0.933, 0.565 )}
colors = [ ( 1., 0.627, 0.478 ), ( 0.565, 0.933, 0.565 ), ( 0.529, 0.808, 0.98 ) ]
colors_dict = {'pulled':( 1., 0.627, 0.478 ), 'semi-pushed':( 0.565, 0.933, 0.565 ), 'fully-pushed':( 0.529, 0.808, 0.98 )}


def cm2inch( x, y ):
    return ( x/2.54, y/2.54 )

def strip( s ): #separates numbers from letters in string
    head = s.strip( '-.0123456789' )
    tail = s[ len( head ): ]
    return head, tail

def get_variables( name ):
    name_root = name.split( '/' )[ -1 ].split( '.' )#Get file name
    if 'txt' in name_root:
        name_root.remove( 'txt' ) #remove ending
    name_root = '.'.join( name_root ) #reform name
    aux = [ strip( s ) for s in name_root.split( '_' ) ]
    #default values if none found
    r0 = 0.01
    m0 = 0.01
    A = 0.0
    B = 0.0
    N = 10000
    for s in aux:
        if s[ 0 ] == 'm':
            m0 = float( s[ 1 ] )
        elif s[ 0 ] == 'A':
            A = float( s[ 1 ] )
        elif s[ 0 ] == 'r':
            r0 = float( s[ 1 ] )
        elif s[ 0 ] == 'B':
            B = float( s[ 1 ] )
        elif s[ 0 ] == 'N':
            N = int( s[ 1 ] )
    return m0, A, r0, B, N


def pick( array, n, offset ):
    ''' Picks n items from array[ offset:len( array ) - offset ], evenly spaced by index'''
    pick_interval = len( array ) - 2*offset
    if n == 1:
        index = int( pick_interval/2 )
    else:
        delta_i = int( pick_interval/( n - 1 ) )
        index = range( offset, offset + pick_interval + 1, delta_i )

    return array[ index ]

def linear_reg( x, y ):
    x = np.asarray( x )
    y = np.asarray( y )
    X = np.vander( x, 2 )
    coeffs, res, rank, sing_vals = np.linalg.lstsq( X, y )
    mx = x.sum()/len( x )
    sx = float( ( ( x - mx )**2 ).sum() )
    if len( x ) > 2:
        r2 = 1. - res/( y.size*np.var( y ) )
    else:
        r2 = 0
    return coeffs, r2


def linear_migration( m0, A, n ):
    return ( m0/2. )*( 1. + A*n )

def pick_color( v, v_F ):
    delta_s = 3./( 2*np.sqrt( 2 ) ) - 1.
    v_semi = ( 1. + delta_s )*v_F
    if v <= v_F:
        color = colors_dict[ 'pulled' ]
    elif v < v_semi:
        color = colors_dict[ 'semi-pushed' ]
    else:
        color = colors_dict[ 'fully-pushed' ]
    return color

def power_law( x, a, b ):
    return a * x ** b

def log_power( x, a, b ):
    return a * ( np.log( x ) ) ** b

def coarse_grain_velocity( v, v_F ):
    delta_s = 3./( 2*np.sqrt( 2 ) ) - 1.
    v_semi = ( 1. + delta_s )*v_F
    v_fully = ( 1. + 2*delta_s )*v_F
    v_out = np.zeros( len( v ) )
    v_out[ np.where( v <= v_F )[ 0 ] ] = v_F
    v_out[ np.where( ( v > v_F )*( v < v_semi ) )[ 0 ] ] = v_semi
    v_out[ np.where( v > v_semi )[ 0 ] ] = v_fully
    return [ v_F, v_semi, v_fully ], v_out

def plot_matrix( x, y, v_array, v_F, coarse=False ):
    if coarse == True:
        A_array = np.zeros( len( x )/10 )
        B_array = np.zeros( len( y )/10 )
        for i in range( len( A_array ) ):
            A_array[ i ] = x[ i*10 ]
        for i in range( len( B_array ) ):
            B_array[ i ] = y[ i*10 ]
    else:
        A_array = x
        B_array = y

    v_values, v_coarse = coarse_grain_velocity( v_array, v_F )
    X, Y = np.meshgrid( A_array, B_array )
    v_matrix = np.zeros( [ len( A_array ), len( B_array ) ] )
    v_plot = np.zeros( [ len( A_array ), len( B_array ) ] )
    for iB in range( len( B_array ) ):
        for iA in range( len( A_array ) ):
            if coarse == True:
                i = ( iB*10 )*len( x ) + iA*10
            else:
                i = iB*len( x ) + iA
            X[ iB ][ iA ] = A_array[ iA ]
            Y[ iB ][ iA ] = B_array[ iB ]
            v_matrix[ iA ][ iB ] = v_array[ i ]
            v_plot[ iA ][ iB ] = v_coarse[ i ]
    return [ X, Y, v_matrix ], v_values, v_plot


def find_fit_endpoint( x_array, y_array, x_fin, epsilon, forward_flag=0 ):
    '''If forward_flag = 1 the search begins at the start of the array
    forward_flag is optional; default value is 0'''
    r_sq_best = 0.

    if forward_flag == 0:
        x_init = x_fin - min( 100, x_fin/2 ) #necessary if nonzero array is small
    else:
        x_init = 0
    x_best = x_init

    while 1 - r_sq_best > epsilon and x_init > 1:
        x_fit = x_array[ x_init:x_fin ]
        y_fit = y_array[ x_init:x_fin ]
        #coeffs, r_sq = linear_reg( x_fit, y_fit )#Do fit
        slope, intercept, r_value, p_value, std_err = stats.linregress( x_fit, y_fit )
        r_sq = r_value ** 2

        if r_sq > r_sq_best:
            x_best = x_init
            r_sq_best = r_sq

        if forward_flag == 0:
            x_init -= 1
        else:
            x_init  += 1

    return x_best, x_array[ x_best:x_fin ], y_array[ x_best:x_fin ]


def fit_heterozygosity( time_arr, het_arr, survival_arr ):
    #calculate fitiing range using sliding window
    #use as final time the point where 5% of simulations have nonzero H
    time_final = sum( het_arr > 0 ) - 1

    # check number of simulations with H != 0; if more than 5%, time_final is last timepoint
    if survival_arr[ time_final ] < 0.05:
        while ( survival_arr[ time_final ] < 0.05 and time_final > 0 ):
            # final point is earliest point with less than 5% of simulations with H != 0
            time_final -= 1

    epsilon = 1E-6 #set threshold for stopping

    time_initial, time_fit, loghet_fit = find_fit_endpoint( time_arr[ :time_final ], np.log( het_arr[ :time_final ] ), time_final, epsilon )

    #coeffs, r_sq = linear_reg( time_fit, np.log( het_fit ) )#Do fit
    slope, intercept, r_value, p_value, std_err = stats.linregress( time_fit, loghet_fit )

    return time_fit, np.exp( loghet_fit ), [ slope, intercept ], r_value ** 2


def read_het( het_data_arr, N, A, B, plot=False, save=True ):
    time_arr = het_data_arr[ 0 ]
    het_arr = het_data_arr[ 1 ]
    survival_arr = het_data_arr[ 2 ]

    time_fit, het_fit, coeffs, r_sq = fit_heterozygosity( time_arr, het_arr, survival_arr )

    neff = 1./abs( coeffs[ 0 ] )

    if ( plot == True or save == True ):
        fig = plt.figure()
        ax = fig.add_subplot( 111 )
        ax.set_title( 'N = ' + str( N ) + ', A = ' + str( A ) + ', B = ' + str( B ) )
        ax.set_xlabel( 'Generations' )
        ax.set_ylabel( 'Log( H )' )

        ax.set_yscale( 'log' )
        ax.plot( time_arr, het_arr, c='0.2', lw=2, alpha=1.0, label='Simulations' )

        fit = np.poly1d( coeffs )
        est = fit( time_fit )
        ax.plot( time_fit, np.exp( est ), c='g', lw=4, alpha=0.8, label='Fit' )
        ax.legend()

        if save == True:
            plt.savefig( output_plots + 'data_plots/heteroplot_N' + str( N ) + '_A' + '{:.3f}'.format( A ) + '_B' + '{:.3f}'.format( B ) + '.pdf' )
    return neff

def contrained_linear( data, param ):
    '''
    returns y = m*x + b
    first row of data should be the constraint:
        [ m, 0 ] costrains slope to m
        [ 0, b ] constrains intercept to b
    '''
    constraint = data[ 0 ]
    x_array = data[ 1: ]
    if constraint[ 0 ] == 0:
        b = constraint[ 1 ]
        return param*x_array + b
    else:
        m = constraint[ 1 ]
        return m*x_array + param

def plot_scaling( Neff_vs_N, coeffs, A, B, figure_fname=None, nu=None ):
    fit = np.poly1d( coeffs )
    est = np.exp( fit( np.log( Neff_vs_N[ :, 0 ] ) ) )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.set_title( 'A = ' + '{:.3f}'.format( A ) + ', B = ' + '{:.3f}'.format( B ) )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.scatter( Neff_vs_N.T[ 0 ], Neff_vs_N.T[ 1 ], c='k' )
    ax.plot( Neff_vs_N.T[ 0 ], est, c='k', label='fit' )

    if nu != None:
        alpha = analytic.exact_exponent( nu )
        const = Neff_vs_N[ len( Neff_vs_N )//2, 1 ]/Neff_vs_N[ len( Neff_vs_N )//2, 0 ]**alpha
        theory = const*Neff_vs_N[ :, 0 ]**alpha
        ax.plot( Neff_vs_N.T[ 0 ], theory, c='r', label='theory' )
        print 'Exponents for A = ', A, ' and B = ', B
        print 'theory: ', alpha, ';  fit: ', coeffs[ 0 ]

    ax.legend( loc='best' )
    if figure_fname != None:
        plt.savefig( figure_fname )


def scaling_analysis( data_directory, plot=False, save_processed_data=True ):
    files_list = glob.glob( data_directory + '/hetero_N*_avg.npy' )

    N_list = []
    A_list = []
    for name in files_list:
        m0, A, r0, B, N = get_variables( name.split( '/' )[ -1 ] )
        N_list.append( N )
        A_list.append( A )
    N_list = np.array( sorted( set( N_list ) ) )
    A_list = np.array( sorted( set( A_list ) ) )

    A_semipushed = 1.5
    comparison_arr = []
    all_data = []
    alpha_array = []
    for A in A_list:
        Neff_arr = []
        files = glob.glob( data_directory + '/hetero_*A' + str( A ) + '_*.npy' )

        for name in files:
            m0, A, r0, B, N = get_variables( name.split( '/' )[ -1 ] )
            #print m0, A, r0, B, N
            hetero_arr = np.load( name )
            Neff = read_het( hetero_arr, N, A, B, plot=True, save=True )
            Neff_arr.append( [ N, Neff, A ] )

            if N == 10000 and ( A == 3 or A == A_semipushed or A == 0 ):
                comparison_arr.append( [ hetero_arr[ 0 ], hetero_arr[ 1 ], hetero_arr[ 2 ] ] )

        Neff_arr = np.array( Neff_arr )
        all_data.append( Neff_arr )

        coeffs, res = linear_reg( np.log( Neff_arr.T[ 0 ][ np.where( Neff_arr.T[ 0 ] >= 10000 ) ] ), np.log( Neff_arr.T[ 1 ][ np.where( Neff_arr.T[ 0 ] >= 10000 ) ] ) )
        if plot == True:
            figure_fname = output_plots + 'data_plots/scaling_A' + '{:.3f}'.format( A ) + '_B' + '{:.3f}'.format( B ) + '.pdf'
            plot_scaling( Neff_arr, coeffs, A, B, figure_fname )
        alpha_array.append( [ A, coeffs[ 0 ] ] )
    alpha_array = np.array( alpha_array )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.scatter( alpha_array.T[ 0 ][ np.where( alpha_array.T[ 0 ] < 50 ) ], alpha_array.T[ 1 ][ np.where( alpha_array.T[ 0 ] < 50 ) ] )
    plt.savefig( output_plots + 'data_plots/linear_scaling.pdf' )

    all_data = np.array( all_data )
    if save_processed_data == True:
        np.save( data_directory + 'neff_data', all_data )
        np.save( data_directory + 'hetero_comparison_data', comparison_arr )

    return all_data


def phaseplot_scaling( data_directory, parameter_file=None, plot=False, plot_het=False, parameter_mapping=None, save_processed_data=True ):
    files_list = sorted( glob.glob( data_directory + '/hetero_N*.npy' ) )
    parameter_dict = {}

    if parameter_mapping != None:
        df_params = pd.read_csv( parameter_mapping[ 0 ] )
        for i, row in df_params.iterrows():
            parameter_dict[ row[ parameter_mapping[ 1 ] ] ] = row[ 'nu' ]

    alpha_array = []
    Neff_array = []
    for fname in files_list:
        m0, A, r0, B, N = get_variables( fname.split( '/' )[ -1 ] )
        hetero_arr = np.load( fname )
        Neff = read_het( hetero_arr, N, A, B, plot=plot_het, save=False )
        Neff_array.append( np.asarray( [ N, Neff, A, B ] ) )
    Neff_array = np.array( Neff_array )

    AB_array = ( Neff_array.T[ [ 2, 3 ] ] ).T
    AB_set = []
    for elem in AB_array:
        if elem.tolist() not in AB_set:
            AB_set.append( elem.tolist() )

    if parameter_file != None:
        parameter_array = np.genfromtxt( parameter_file, delimiter=',' )

    output_array = []
    for elem in AB_set:
        # Get data from Neff_array
        [A, B] = elem
        Neff_subset = Neff_array[ np.where( Neff_array.T[2] == A )[ 0 ] ]
        Neff_subset = Neff_subset[ np.where( Neff_subset.T[3] == B )[ 0 ] ]

        if parameter_file != None:
            if ( A in parameter_array[ :, 0 ] ) and ( B in parameter_array[ :, 1 ] ):
                output_array.append( Neff_subset )
        else:
            output_array.append( Neff_subset )

        # Sort by N
        Neff_subset = Neff_subset[ Neff_subset[ :, 0 ].argsort() ]
        slope, intercept, r_value, p_value, std_err = stats.linregress( np.log( Neff_subset[ :, 0 ] ), np.log( Neff_subset[ :, 1 ] ) )
        if A in parameter_dict:
            nu = parameter_dict[ A ]
        else:
            nu = None

        if plot == True:
            figure_fname = output_plots + 'data_plots/linear_A' + '{:.3f}'.format( A ) + '{:.3f}'.format( B ) + '.pdf'
            if nu != 1:
                plot_scaling( Neff_subset, [ slope, intercept ], A, B, figure_fname, nu=nu )
            else:
                plot_scaling( Neff_subset, [ slope, intercept ], A, B, figure_fname )
        alpha_array.append( [ A, B, slope ] )

    alpha_array = np.array( alpha_array )
    output_array = np.array( output_array )

    # Plot alpha vs A
    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.scatter( alpha_array.T[ 0 ], alpha_array.T[ 1 ] )
    plt.savefig( output_plots + 'data_plots/phaseplot_vs_A.pdf' )

    if save_processed_data == True:
        np.save( data_directory + 'alpha_phaseplot', alpha_array )
        np.save( data_directory + 'neff_data_phaseplot', output_array )

    return output_array


def check_scaling( data_directory, ab_array=[ [ 1.371, 0.546 ] ], rm_array=[ [ 0.01, 0.01 ] ], plot=True, save_processed_data=False ):

    for ab in ab_array:
        [ A, B ] = ab

        for rm in rm_array:
            alpha_array = []
            Neff_array = []
            [ r0, m0 ] = rm
            files_list = sorted( glob.glob( data_directory + '/hetero_N*_r' + str( r0 ) + '_m' + str( m0 ) + '_*.npy' ) )

            for fname in files_list:
                m0, A, r0, B, N = get_variables( fname.split( '/' )[ -1 ] )
                hetero_arr = np.load( fname )
                Neff = read_het( hetero_arr, N, A, B, plot=False, save=False )
                Neff_array.append( np.asarray( [ N, Neff, A, B ] ) )

            Neff_array = np.array( Neff_array )

            AB_array = ( Neff_array.T[ [ 2, 3 ] ] ).T
            AB_set = []
            for elem in AB_array:
                if elem.tolist() not in AB_set:
                    AB_set.append( elem.tolist() )

            output_array = []
            for elem in AB_set:
                # Get data from Neff_array
                [A, B] = elem
                Neff_subset = Neff_array[ np.where( Neff_array.T[2] == A )[ 0 ] ]
                Neff_subset = Neff_subset[ np.where( Neff_subset.T[3] == B )[ 0 ] ]
                output_array.append( Neff_subset )

                # Sort by N
                Neff_subset = Neff_subset[ Neff_subset[ :, 0 ].argsort() ]
                slope, intercept, r_value, p_value, std_err = stats.linregress( np.log( Neff_subset[ :, 0 ] ), np.log( Neff_subset[ :, 1 ] ) )
                #fit = np.poly1d( [ slope, intercept ] )
                #est = np.exp( fit( np.log( Neff_subset.T[ 0 ] ) ) )

                if plot == True:
                    figure_fname = output_plots + 'data_plots/linear_frontwidth_A' + '{:.3f}'.format( A ) + '{:.3f}'.format( B ) + '.pdf'
                    plot_scaling( Neff_subset, [ slope, intercept ], A, B, figure_fname )
                alpha_array.append( [ A, B, slope ] )
                #print A, B, 2*A + B, slope

            alpha_array = np.array( alpha_array )
            output_array = np.array( output_array )

            # Plot alpha vs A
            fig = plt.figure()
            ax = fig.add_subplot( 111 )
            ax.scatter( alpha_array.T[ 0 ], alpha_array.T[ 1 ] )
            plt.savefig( output_plots + 'data_plots/front_width_vs_A.pdf' )

            '''
            if save_processed_data == True:
                np.save( data_directory + 'alpha_phaseplot', alpha_array )
                np.save( data_directory + 'neff_data_phaseplot', output_array )
            '''

    return output_array


def profile( data_directory ):
    files_list = glob.glob( data_directory + '/profile_*.npy' )

    N_list = []
    m0_list = []
    for name in files_list:
        m0, A, r0, B, N = get_variables( name.split( '/' )[ -1 ] )
        N_list.append( N )
        m0_list.append( m0 )
    N_list = np.array( sorted( set( N_list ) ) )
    m0_list = np.array( sorted( set( m0_list ) ) )

    for m0 in m0_list:
        Neff_arr = []
        files = glob.glob( data_directory + '/profile_*m0-' + str( m0 ) + '_*.npy' )

        for name in files:
            m0, A, r0, B, N = get_variables( name.split( '/' )[ -1 ] )
            if m0 != 0:
                A = ( 0.25 - m0 )/m0
            else:
                A = 100

            prof_data = np.load( name )
            x = prof_data.T[ 0 ]
            prof_arr = ( prof_data.T[ 1 ] + prof_data.T[ 2 ] )/float( N )

            #fig = plt.figure()
            #ax = fig.add_subplot( 121 )
            #ax.plot( x, prof_arr )
            #ax = fig.add_subplot( 122 )
            #ax.set_yscale( 'log' )
            #ax.plot( x[ np.nonzero( prof_arr )[ 0 ] ], prof_arr[ np.nonzero( prof_arr )[ 0 ] ] )

            #plt.savefig( output_plots + 'plots/' + data_directory + '_profile_N' + str( N ) + '_mr' + '{:.3f}'.format( A ) + '.pdf' )


def velocity( data_directory ):
    files_list = glob.glob( data_directory + '/velocity_*.npy' )

    N_list = []
    m0_list = []
    for name in files_list:
        m0, A, r0, B, N = get_variables( name.split( '/' )[ -1 ] )
        N_list.append( N )
        m0_list.append( m0 )
    N_list = np.array( sorted( set( N_list ) ) )
    m0_list = np.array( sorted( set( m0_list ) ) )

    alpha_array = []
    for m0 in m0_list:
        Neff_arr = []
        files = glob.glob( data_directory + '/velocity_*m0-' + str( m0 ) + '_*.npy' )

        for name in files:
            m0, A, r0, B, N = get_variables( name.split( '/' )[ -1 ] )
            if m0 != 0:
                A = ( 0.25 - m0 )/m0
            else:
                A = 100
            vel_arr = np.load( name )
            coeffs, res = linear_reg( vel_arr.T[ 0 ][ len( vel_arr )/2: ], vel_arr.T[ 1 ][ len( vel_arr )/2: ] )
            Neff_arr.append( [ N, coeffs[ 0 ], 2*np.sqrt( m0*0.01/2. ) ] )
        Neff_arr = np.array( Neff_arr )

        coeffs, res = linear_reg( Neff_arr.T[ 0 ][ len( Neff_arr )/2: ], Neff_arr.T[ 1 ][ len( Neff_arr )/2: ] )
        fit = np.poly1d( coeffs )
        est = fit( Neff_arr.T[ 0 ] )

        fig = plt.figure()
        ax = fig.add_subplot( 111 )
        ax.set_title( 'm1/m0 = ' + '{:.3f}'.format( A ) )
        ax.scatter( Neff_arr.T[ 0 ], Neff_arr.T[ 1 ] )
        ax.plot( Neff_arr.T[ 0 ], est )
        plt.savefig( output_plots + 'data_plots/' + data_directory + '_velocit_mr' + '{:.3f}'.format( A ) + '.pdf' )

        alpha_array.append( [ A, coeffs[ 0 ] ] )
    alpha_array = np.array( alpha_array )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.scatter( alpha_array.T[ 0 ][ np.where( alpha_array.T[ 0 ] < 50 ) ], alpha_array.T[ 1 ][ np.where( alpha_array.T[ 0 ] < 50 ) ] )
    plt.savefig( output_plots + 'data_plots/' + data_directory + '_scaling.pdf' )


def plot_diversity( data_arr, m_pushed, m_semi, m_pulled ):
    font = {'family' : 'sans-serif', 'serif' : 'Helvetica Neue', 'weight' : 'bold', 'size' : 12}
    matplotlib.rc( 'font', **font )

    het_data = np.load( input_path + 'hetero_comparison_data.npy' )
    axis_fontsize = 12

    fig = plt.figure( figsize=( cm2inch( 17.8 ), cm2inch( 9.0 ) ) )
    ax = fig.add_subplot( 121 )
    ax.set_yscale( 'log' )
    ax.ticklabel_format( style='sci', scilimit=( -2,2 ), axis='x' )
    ax.set_xlabel( 'time, t', fontweight='bold', fontsize=axis_fontsize )
    ax.set_ylabel( 'heterozygosity, H', fontweight='bold', fontsize=axis_fontsize )
    ax.text( -15000, 1.5E1, 'A', fontweight='bold', fontsize=14 )
    ax.text( 70000, 1.5E-1, '$H \sim e^{-t/T_c}$', fontsize=14, color='k' )
    ax.set_xticks( [ 0, 40000, 80000 ] )
    ax.set_yticks( [ 1E-4, 1E-2, 1 ] )
    ax.set_xlim( [ 0, 120000 ] )
    ax.set_ylim( [ 1E-4, 1E1 ] )

    number_of_points = 50
    reduced_het_indices = ( 1000/number_of_points )*np.arange( number_of_points )

    for index, elem in enumerate( het_data ):
        x_het = elem[ 0 ]
        het = elem[ 1 ]
        surv = elem[ 2 ]

        x_plot = np.array( [ x_het[ i ] for i in reduced_het_indices ] )
        het_plot = np.array( [ het[ i ] for i in reduced_het_indices ] )

        het_plot = het_plot[ np.where( x_plot < 100000 )[ 0 ] ]
        x_plot = x_plot[ np.where( x_plot < 100000 )[ 0 ] ]

        xh, het_fit, coeffs, r_sq = fit_heterozygosity( x_het, het, surv )
        fit = np.poly1d( coeffs )
        x_fit = [ x_het[ 50 ], 1.1*x_het[ -1 ] ] #choose range for plotting fit
        est = np.exp( fit( x_fit ) )

        #Plot results
        if index==0:
            clr = 'b'
            lbl = 'pushed'
        else:
            clr = 'r'
            lbl = 'pulled'
        ax.scatter( x_plot, het_plot, s=30, edgecolor=clr, facecolor='none', lw=1, label=lbl )
        ax.plot( x_fit, est, c=clr, lw=1 )

    legend_properties={'weight':'bold', 'size':10}
    ax.legend( loc='best', prop=legend_properties, scatterpoints=1 )


    ax = fig.add_subplot( 122 )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    #ax.ticklabel_format( style='sci', scilimit=( -2,2 ), axis='x' )
    ax.set_xlabel( 'population size, N', fontweight='bold', fontsize=axis_fontsize )
    ax.set_ylabel( 'coalescence time, $\mathbf{T_c}$', fontweight='bold', fontsize=axis_fontsize )
    ax.text( 2E2, 2.5E7, 'B', fontweight='bold', fontsize=14 )
    ax.set_xticks( [ 1E3, 1E5, 1E7 ] )
    ax.set_yticks( [ 1E3, 1E5, 1E7 ] )
    ax.set_xlim( [ 7E2, 1E7 ] )
    ax.set_ylim( [ 7E2, 2E7 ] )

    for index in range( n_examples ):
        n_arr = data_arr[ m_pushed[ index ] ].T[ 0 ]
        tc_arr = data_arr[ m_pushed[ index ] ].T[ 1 ]
        if index == 0:
            ax.scatter( n_arr, tc_arr, edgecolor='none', facecolor='blue', s=30, label='fully-pushed' )
        else:
            ax.scatter( n_arr, tc_arr, edgecolor='none', facecolor='blue', s=30 )

        n_arr = data_arr[ m_semi[ index ] ].T[ 0 ]
        tc_arr = data_arr[ m_semi[ index ] ].T[ 1 ]
        if index == 0:
            ax.scatter( n_arr, tc_arr, edgecolor='none', facecolor='green', s=30, label='semi-pushed' )
        else:
            ax.scatter( n_arr, tc_arr, edgecolor='none', facecolor='green', s=30 )

        n_arr = data_arr[ m_pulled[ index ] ].T[ 0 ]
        tc_arr = data_arr[ m_pulled[ index ] ].T[ 1 ]
        if index == 0:
            ax.scatter( n_arr, tc_arr, edgecolor='none', facecolor='red', s=30, label='pulled' )
        else:
            ax.scatter( n_arr, tc_arr, edgecolor='none', facecolor='red', s=30 )

    legend_properties={'weight':'bold', 'size':10}
    ax.legend( loc='upper left', prop=legend_properties, scatterpoints=1 )

    plt.tight_layout( pad=2.0 )
    plt.savefig( output_plots + 'data_plots/diversity_loss.pdf' )


def plot_models( full_data_arr, m_arr ):
    font = {'family' : 'sans-serif', 'serif' : 'Helvetica Neue', 'weight' : 'bold', 'size' : 12}
    matplotlib.rc( 'font', **font )

    axis_fontsize = 12
    n_examples = len( m_arr )

    fig = plt.figure( figsize=( cm2inch( 11.4 ), cm2inch( 9.0 ) ) )
    ax = fig.add_subplot( 111 )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    #ax.ticklabel_format( style='sci', scilimit=( -2,2 ), axis='x' )
    ax.set_xlabel( 'population size, N', fontweight='bold', fontsize=axis_fontsize )
    ax.set_ylabel( 'coalescence time, $\mathbf{T_c}$', fontweight='bold', fontsize=axis_fontsize )
    #ax.text( 2E2, 2.5E7, 'B', fontweight='bold', fontsize=14 )
    ax.set_xticks( [ 1E3, 1E5, 1E7 ] )
    ax.set_yticks( [ 1E3, 1E5, 1E7 ] )
    ax.set_xlim( [ 7E2, 1E7 ] )
    ax.set_ylim( [ 7E2, 2E7 ] )

    markers = [ 'o', 's', '^' ]
    labels = [ '$m_0 + m_1 \\rho$', '$m_0 + m_2 \\rho^2$', '$m_0 + m_{1/2} \sqrt{\\rho}$' ]
    colors = [ 'blue', 'green', 'red' ]
    c_model = 0
    for data_arr in full_data_arr:
        for index in range( n_examples ):
            n_arr = data_arr[ m_arr[ index ] ].T[ 0 ]
            tc_arr = data_arr[ m_arr[ index ] ].T[ 1 ]
            if index == 0:
                ax.scatter( n_arr, tc_arr, marker=markers[ c_model ], edgecolor=colors[ index ], facecolor='none', s=30, label=labels[ c_model ] )
            else:
                ax.scatter( n_arr, tc_arr, marker=markers[ c_model ], edgecolor=colors[ index ], facecolor='none', s=30 )
        c_model  += 1


    legend_properties={'weight':'bold', 'size':10}
    ax.legend( loc='upper left', prop=legend_properties, scatterpoints=1 )

    plt.tight_layout( pad=2.0 )
    plt.savefig( output_plots + 'data_plots/migration_models.pdf' )


def talk_plot( data_arr, m_pushed, m_semi, m_pulled ):
    font = {'family' : 'sans-serif', 'serif' : 'Helvetica Neue', 'weight' : 'bold', 'size' : 12}
    matplotlib.rc( 'font', **font )

    het_data = np.load( 'hetero_comparison_data.npy' )
    axis_fontsize = 12

    fig = plt.figure( figsize=( cm2inch( 11.4 ), cm2inch( 9.0 ) ) )
    ax = fig.add_subplot( 111 )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    #ax.ticklabel_format( style='sci', scilimit=( -2,2 ), axis='x' )
    ax.set_xlabel( 'population size, N', fontweight='bold', fontsize=axis_fontsize )
    ax.set_ylabel( 'rate of diversity loss, $\mathbf{\Lambda}$', fontweight='bold', fontsize=axis_fontsize )
    ax.set_xticks( [ 1E4, 1E5, 1E6 ] )
    ax.set_yticks( [ 1E-3, 1E-5, 1E-7 ] )
    ax.set_xlim( [ 5E3, 1E7 ] )
    ax.set_ylim( [ 8E-8, 2E-3 ] )

    for index in range( n_examples ):
        n_arr = data_arr[ m_pushed[ index ] ].T[ 0 ]
        tc_arr = data_arr[ m_pushed[ index ] ].T[ 1 ]
        n_fit = [ 0.5*min( n_arr ), 2*max( n_arr ) ]
        fit = ( n_arr[ -1 ]/tc_arr[ -1 ] )/n_fit
        if index == 0:
            ax.scatter( n_arr, 1./tc_arr, edgecolor='blue', facecolor='none', s=50, label='fully-pushed' )
            ax.plot( n_fit, fit, lw=1, c='blue' )
        else:
            ax.scatter( n_arr, 1./tc_arr, edgecolor='blue', facecolor='none', s=50 )
            ax.plot( n_fit, fit, lw=1, c='blue' )

        n_arr = np.array( sorted( data_arr[ m_semi[ index ] ].T[ 0 ] ) )
        tc_arr = np.array( sorted( data_arr[ m_semi[ index ] ].T[ 1 ] ) )
        n_fit = [ 0.5*min( n_arr ), 2*max( n_arr ) ]
        coeffs, res = linear_reg( np.log( n_arr[ 1: ] ), -np.log( tc_arr[ 1: ] ) )
        fit = np.poly1d( coeffs )
        est = np.exp( fit( np.log( n_fit ) ) )
        if index == 0:
            ax.scatter( n_arr, 1./tc_arr, edgecolor='green', facecolor='none', s=50, label='semi-pushed' )
            ax.plot( n_fit, est, lw=1, c='green' )
        else:
            ax.scatter( n_arr, 1./tc_arr, edgecolor='green', facecolor='none', s=50 )
            ax.plot( n_fit, est, lw=1, c='green' )

        n_arr = data_arr[ m_pulled[ index ] ].T[ 0 ]
        tc_arr = data_arr[ m_pulled[ index ] ].T[ 1 ]
        n_fit = [ 0.5*min( n_arr ), n_arr[ 2 ], n_arr[ 4 ], 2*max( n_arr ) ]
        fit = np.log( sorted( n_fit ) )**( -3 )*( np.log( n_arr[ -1 ] )**3/tc_arr[ -1 ] )
        if index == 0:
            ax.scatter( n_arr, 1./tc_arr, edgecolor='red', facecolor='none', s=50, label='pulled' )
            ax.plot( sorted( n_fit ), fit, lw=1, c='red' )
        else:
            ax.scatter( n_arr, 1./tc_arr, edgecolor='red', facecolor='none', s=50 )
            ax.plot( sorted( n_fit ), fit, lw=1, c='red' )

        ax.text( 2.5E6, 8E-5, '$log^{-3}N$', color='red' )
        ax.text( 3.0E6, 8E-6, '$N^{\\alpha}$', color='green' )
        ax.text( 3.0E6, 2E-7, '$N^{-1}$', color='blue' )

    legend_properties={'weight':'bold', 'size':10}
    #ax.legend( loc='upper right', prop=legend_properties, scatterpoints=1 )

    plt.tight_layout( pad=2.0 )
    plt.savefig( output_plots + 'data_plots/talk_density_dependence.pdf' )


def check_simulation_time():
    parameters = np.loadtxt( "2d_parameters.csv", delimiter="," )
    pushed_boundary = 2.0
    fullypushed_boundary = 4.0
    m0 = 0.01
    r0 = 0.01
    for param in parameters:
        [ A, B ] = param
        if 2*A + B >= pushed_boundary and 2*A + B < fullypushed_boundary:
            exponent = ( 2*A + B  - pushed_boundary )/( fullypushed_boundary - pushed_boundary )
            #print A, B, exponent
            for N in [10**4, 10**5, 10**6]:
                t_semipushed = np.sqrt( m0/( 2*r0 ) )*N**exponent
                #print N, 20*t_semipushed, 200*np.sqrt( N )
            #print '\n'


def check_parameters( A, B ):
    if 2*A + B <= 2.0:
        return 'pulled'
    elif 2*A + B < 4.0:
        return 'semi-pushed'
    else:
        return 'fully-pushed'


def measure_velocity( data_path, parameter_file=None, plot=True ):
    if parameter_file != None:
        parameter_array = np.genfromtxt( parameter_file, delimiter=',' )
        v_pulled = []
        v_semi = []
        v_fully = []
        for [ A, B ] in parameter_array:
            velocity_files = glob.glob( data_path + '*_A' + str( A ) + '_B' + str( B ) + '_*.npy' )

            # Get largest N
            _, _, _, _, Nmax = get_variables( velocity_files[ 0 ] )
            for fname in velocity_files:
                m0, A, r0, B, N = get_variables( fname )
                Nmax = max( N, Nmax )

            velocity_data = np.load( data_path + 'velocity_N' + str( Nmax ) + '_r' + str( r0 ) + '_m' + str( m0 ) + '_A' + str( A ) + '_B' + str( B ) + '_avg.npy' )

            # Exclude first 10% of simulations
            t_array = velocity_data[ 100:, 0 ]
            n_array = velocity_data[ 100:, 1 ]

            dt_array = np.array( [ t_array[ i ] - t_array[ i + 1 ] for i in range( len( t_array ) - 1 ) ] )
            dn_array = np.array( [ n_array[ i ] - n_array[ i + 1 ] for i in range( len( n_array ) - 1 ) ] )
            v_mean = np.mean( dn_array/dt_array )
            v_slope = ( n_array[ -1 ] - n_array[ 0 ] )/( t_array[ -1 ] - t_array[ 0 ] )
            v_regression = stats.linregress( t_array, n_array )[ 0 ]

            if check_parameters( A, B ) == 'pulled':
                v_pulled.append( v_slope )
            elif check_parameters( A, B ) == 'semi-pushed':
                v_semi.append( v_slope )
            elif check_parameters( A, B ) == 'fully-pushed':
                v_fully.append( v_slope )

        #print v_pulled
        #print v_semi
        #print v_fully

        v_pulled_avg = np.mean( v_pulled )
        v_semi_avg = np.mean( v_semi )
        v_fully_avg = np.mean( v_fully )

        #print 'Pulled: nu = ', v_pulled_avg/v_pulled_avg, ', v = ', v_pulled_avg
        #print 'Semi-pushed: nu = ', v_semi_avg/v_pulled_avg, ', v = ', v_semi_avg
        #print 'Fully-pushed: nu = ', v_fully_avg/v_pulled_avg, ', v = ', v_fully_avg, '\n'

        #print 'Pulled exponent: ', analytic.exact_exponent( v_pulled_avg/v_pulled_avg )
        #print 'Semi-pushed exponent: ', analytic.exact_exponent( v_semi_avg/v_pulled_avg )
        #print 'Fully-pushed exponent: ', analytic.exact_exponent( v_fully_avg/v_pulled_avg )


def find_parameters( nu_array, epsilon=1E-4, save_file=None ):
    phase_plot_file = '../data/data_2d_linear-migration_cooperative-growth.txt'

    # Define constants
    r0 = 0.001
    m0 = 0.05
    v_F = 2.*np.sqrt( r0*m0/2. )

    # Load data
    data_array = np.loadtxt( phase_plot_file, delimiter=',' )

    x = np.sort( np.unique( data_array.T[ 0 ] ) )
    y = np.sort( np.unique( data_array.T[ 1 ] ) )
    v_array = data_array.T[ 2 ]
    [ x_matrix, y_matrix, v_matrix ], v_values, v_plot = plot_matrix( x, y, v_array, v_F )


    # Choose parameters
    params_array = [ data_array[ np.where( abs( v_array - nu*v_F )/( nu*v_F ) < epsilon )[ 0 ] ] for nu in nu_array ]
    one2one_array = np.array( [ pick( params, 1, int( len( params )*0.25 ) ) for params in params_array ] )

    if save_file != None:
        np.savetxt( save_file, one2one_array[ :, [ 0, 1 ] ], fmt='%.3f', delimiter=',' )


def get_fig4_parameters( nu_array, epsilon=1E-4, save_file=None ):
    phase_plot_file = '../data/data_2d_linear-migration_cooperative-growth.txt'

    # Define constants
    r0 = 0.001
    m0 = 0.05
    v_F = 2.*np.sqrt( r0*m0/2. )

    # Load data
    data_array = np.loadtxt( phase_plot_file, delimiter=',' )

    x = np.sort( np.unique( data_array.T[ 0 ] ) )
    y = np.sort( np.unique( data_array.T[ 1 ] ) )
    v_array = data_array.T[ 2 ]
    [ x_matrix, y_matrix, v_matrix ], v_values, v_plot = plot_matrix( x, y, v_array, v_F )

    v_dispersal = v_matrix[ np.where( y_matrix == 0 ) ]
    entries = len( v_dispersal )
    index_array = [ np.argmin( abs( v_dispersal - v_F*nu ) ) for nu in nu_array ]

    A_param = x_matrix[ np.where( y_matrix == 0 ) ][ index_array ]
    v_param = v_dispersal[ index_array ]

    if save_file == None:
        return A_array

    with open( save_file, 'w' ) as fout:
        fout.write( 'nu\tv\tA linear\tA sq. root\tA quadratic\n' )

        for i, nu in enumerate( nu_array ):
            index = index_array[ i ]
            if index == entries - 1:
                A = ''
                v = ''
            else:
                v = str( v_param[ i ] )
                A = str( A_param[ i ] )
            #fout.write( str( nu ) + '\t' + v + '\t' + A + '\t' + '' + '\t' + '' + '\n' )
            fout.write( "{0:.3f}\t{1}\t{2}\t{3}\t{4}\n".format( nu, v, A, '', '' ) )


def fraction_distribution( data_dir, N, A=0.25, B=0.0, p_target=None, epsilon_p=0.01, n_plots='all', bins=30, save_plots=True ):
    file_list = sorted( glob.glob( data_dir + '/fraction_*_N' +str( N ) + '_A' + str( A ) + '_B' + str( B ) + '_points=auto.csv' ) )
    wave_dict = { 0:'pulled', 1:'semi-pushed', 2:'fully-pushed' }

    for fname in file_list:
        _, B, _, A, N = get_variables( fname )
        wave_type = wave_dict[ analytic.linear_cooperativity( A, B ) ]
        identifier = 'N' + str( N ) + '_A' + str( A ) + '_B' + str( B )

        fraction_array = np.loadtxt( fname, delimiter=',' )
        ( time_samples, runs ) = fraction_array.shape
        if ( n_plots == 'all' ) or ( p_target != None ):
            time_index = range( time_samples )
        else:
            time_index = np.linspace( 0, time_samples - 1, n_plots, dtype=int )

        for i in time_index:
            f_survived =  fraction_array[ i ][ np.where( ( fraction_array[ i ] != 0.0 )*( fraction_array[ i ] != 1.0 ) )[ 0 ] ]
            p_survived = len( f_survived )/float( runs )

            # If given p_target, skip plots that don't satisfy condition
            if p_target != None:
                if abs( p_survived - p_target ) > epsilon_p:
                    continue

            fig = plt.figure()
            ax = fig.add_subplot( 111 )
            ax.set_title( 'N = ' + str( N ) + ', ' + wave_type + ', $p_{surv} = $' + str( p_survived )  )
            ax.set_xlabel( 'allele fraction, f' )
            ax.set_ylabel( 'probability density' )
            ax.hist( f_survived, bins=bins, density=True, color='grey', label='simulations' )

            if ( wave_type == 'pulled' ) or ( wave_type == 'semi-pushed' ):
                ax.set_yscale( 'log' )
                hist = np.histogram( f_survived, bins=100, density=True )
                p_density = hist[ 0 ]
                c = p_density[ int( 2*len( p_density )/5 ):int( 3*len( p_density )/5 ) ].mean()

                if c == 0.0:
                    c = min( p_density[ np.where( p_density != 0.0 )[ 0 ] ] )

                f_plot = np.linspace( 1E-2, 1.0 - 1E-2, 1000 )
                p_plot = ( c/4 )/( f_plot*( 1.0 - f_plot ) )
                ax.plot( f_plot, p_plot, lw=2, c='r', label='$\propto \\frac{1}{f(1-f)}$' )

                if wave_type == 'semi-pushed':
                    alpha=0.5
                    p_sqroot_plot = ( c/4**alpha )/( ( f_plot*( 1.0 - f_plot ) )**alpha )
                    ax.plot( f_plot, p_sqroot_plot, lw=2, c='g', ls=':', label='$\propto \\frac{1}{ [ f(1-f) ]^{0.5}}$' )

                    alpha=0.8
                    p_sqroot_plot = ( c/4**alpha )/( ( f_plot*( 1.0 - f_plot ) )**alpha )
                    ax.plot( f_plot, p_sqroot_plot, lw=2, c='g', ls='--', label='$\propto \\frac{1}{ [ f(1-f) ]^{0.8}}$' )

                #if wave_type == 'pulled':
                #    p_corrected_plot = ( c/4 )/( -f_plot*np.log( 1.0 - f_plot )*( f_plot - 1.0 )*np.log( 1.0 - f_plot ) )
                #    ax.plot( f_plot, p_plot, lw=2, ls='--', c='purple', label='$\propto \\frac{1}{f(1 - f) \log f\log(1 - f)}$' )

                ax.legend( loc='upper center' )

            if save_plots == True:
                plt.savefig( output_plots + 'fraction_distributions/fraction_' + identifier + '_p' + str( p_survived ) + '_bins' + str( bins ) + '.pdf' )
            plt.close()


def fraction_distribution_all( data_dir, save_plots=True ):
    file_list = sorted( glob.glob( data_dir + '/fraction_*_points=auto.csv' ) )
    wave_dict = { 0:'pulled', 1:'semi-pushed', 2:'fully-pushed' }

    for fname in file_list:
        _, B, _, A, N = get_variables( fname )
        #print N, A, B
        wave_type = wave_dict[ analytic.linear_cooperativity( A, B ) ]

        fraction_array = np.loadtxt( fname, delimiter=',' )
        ( time_samples, runs ) = fraction_array.shape
        time_index = range( time_samples )
        #time_index = [ 0 ]

        fig = plt.figure()
        ax = fig.add_subplot( 111 )
        ax.set_title( 'N = ' + str( N ) + ', ' + wave_type  )
        ax.set_xlabel( 'allele fraction, f' )
        ax.set_ylabel( 'probability density' )

        delta_survival = 1.0
        i_target = 0
        p_target = 0.5
        for i in time_index:
            f_survived =  fraction_array[ i ][ np.where( ( fraction_array[ i ] != 0.0 )*( fraction_array[ i ] != 1.0 ) )[ 0 ] ]
            p_survived = len( f_survived )/float( runs )
            delta = abs( p_survived - p_target )

            if delta < delta_survival:
                i_target = i
                delta_survival = delta

        f_survived =  fraction_array[ i_target ][ np.where( ( fraction_array[ i_target ] != 0.0 )*( fraction_array[ i_target ] != 1.0 ) )[ 0 ] ]
        #print 'fraction survived = ', len( f_survived )/float( runs )
        ax.hist( f_survived, bins=50, density=True, label='time=' + str( i_target ) )
        ax.legend( loc='best' )

        #print '\n'


def parse_profile_data( profile_data, index ):
    data_t = profile_data[ index ]
    t = data_t[ 0 ]
    x_array = data_t[ 1:301 ]
    n1_array = data_t[ 301:601 ]
    n2_array = data_t[ 601: ]
    N = max( n1_array + n2_array )
    rho1_array = n1_array/N
    rho2_array = n2_array/N

    return t, N, x_array, rho1_array, rho2_array


def read_fraction_distribution( fname, n_cutoff ):
    data_file = np.loadtxt( fname, delimiter=',' )
    l_box = ( data_file.shape[ 1 ] - 1 ) / 2
    t_array = data_file[ :, 0 ]
    n1_array = np.array( data_file[ :, [ i for i in range( 1, 2 * l_box, 2 ) ] ] )
    n2_array = np.array( data_file[ :, [ i for i in range( 2, 2 * l_box + 1, 2 ) ] ] )
    n_array = n1_array + n2_array
    f_array = n1_array
    f_array[ np.where( n_array != 0 ) ] = f_array[ np.where( n_array != 0 ) ] / n_array[ np.where( n_array != 0 ) ]

    f_mean = []
    f_std = []
    for i, t in enumerate( t_array ):
        i_max = np.arange( l_box )[ np.where( n_array[ i ] > n_cutoff )[ 0 ] ][ -1 ] + 1
        f_m = ( f_array[ i ] * n_array[ i ] ).sum() / n_array[ i ].sum()
        f_mean.append( f_m )
        f_s = np.sqrt( ( ( f_array[ i ][ :i_max ] - f_mean[ -1 ] ) ** 2 ).mean() )
        f_std.append( f_s )
    f_mean = np.array( f_mean )
    f_std = np.array( f_std )

    return t_array, n_array, f_array, f_mean, f_std


def group_velocities():
    # Group velocity data
    file_list = glob.glob( input_path + 'data_cooperation/r' + str( r_plot ) + '_m' + str( m_plot ) + '/velocity_N' + str( N_plot ) + '_*B0.0_*.npy' )
    velocity_data = []
    for fname in file_list:
        _, A, _, _, _ = analysis.get_variables( fname )
        v_array = np.load( fname )

        if v_measurement == 'ratio':
            i_index = int( 0.1*len( v_array ) )
            f_index = -1
            v = ( v_array[ f_index, 1 ] - v_array[ i_index, 1 ] )/( v_array[ f_index, 0 ] - v_array[ i_index, 0 ] )
            velocity_data.append( [ A, v ] )
    velocity_data = np.array( velocity_data )
    velocity_data = velocity_data[ np.argsort( velocity_data[ :, 0 ] ) ]

    return velocity_data


if __name__=='__main__':
    #linear_neff = scaling_analysis( 'linear_data', 0, 1 )
    #linear_neff = scaling_analysis( 'quadratic_data', 0, 1 )
    #linear_neff = scaling_analysis( 'sqroot_data', 0, 1 )

    linear_neff = np.load( input_path + 'neff_linear_data.npy' )
    quadratic_neff = np.load( input_path + 'neff_quadratic_data.npy' )
    sqroot_neff = np.load( input_path + 'neff_sqroot_data.npy' )

    data_arr = linear_neff
    n_examples = 2
    A = linear_neff.T[ 2 ][ 0 ]
    m_pushed = [ 0, 2 ]
    m_semi = [ 4, 7 ]
    m_pulled = [ -1, -2 ]
    m_comparison = [ 0, 4, -1 ]

    #plot_diversity( linear_neff, m_pushed, m_semi, m_pulled )
    #plot_models( [ linear_neff, quadratic_neff, sqroot_neff ], m_comparison )
    #plot_diversity( quadratic_neff, m_pushed, m_semi, m_pulled )
    #plot_diversity( sqroot_neff, m_pushed, m_semi, m_pulled )

    #talk_plot( linear_neff, m_pushed, m_semi, m_pulled )
    #scaling_analysis( input_path + 'phaseplot_data/run3_test/', plot=True, save_processed_data=True )
    #scaling_analysis( input_path + 'data_front_width/', plot=True, save_processed_data=False )
    #phaseplot_scaling( input_path + 'data_cooperation/r0.001_m0.05/', plot=True, save_processed_data=True )
    model = 'sqrt'
    phaseplot_scaling( input_path + 'scaling_data/' + model + '_model/', plot=True, plot_het=False, parameter_mapping=[ '../run_parameters/fig4_parameters.csv', 'm1/m0_' + model ], save_processed_data=True )
    #check_scaling( input_path + 'data_front_width/', rm_array=[ [ 0.01, 0.01 ], [ 0.001, 0.05 ], [ 0.0001, 0.05 ] ], plot=True, save_processed_data=False )
    #print( analytic.exact_exponent( 1.05 ) )

    #measure_velocity( input_path + 'data_phaseplot/', '2d_parameters.csv' )

    #find_parameters( [ 1.01, 1.02, 1.03, 1.04, 1.05, 1.055, 1.06, 1.08, 1.1, 1.15, 1.2, 1.25, 1.3, 1.5 ], save_file='2d_full_scaling_parameters.csv' )
    #delta_nu = 3/( 2*np.sqrt( 2 ) ) - 1
    #get_fig4_parameters( [ 1 + delta_nu/2, 1 + 2*delta_nu/3, 1 + 4*delta_nu/5, 1 + delta_nu, 1.1, 1.2, 1.3, 1.5 ], save_file='scaling_parameters.csv' )


    #profile( 'linear_data' )
    #profile( 'quadratic_data' )
    #check_simulation_time()
    #A = 1.44
    #B = 2.922
    A = 0.693
    B = 1.23
    #fraction_distribution( input_path + '/data_fraction_distribution', N=10000, A=A, B=B, p_target=0.4 )
    #fraction_distribution( input_path + '/data_fraction_distribution', N=100000, A=A, B=B, p_target=0.4 )
    #fraction_distribution( input_path + '/data_fraction_distribution', N=1000000, A=A, B=B, p_target=0.4 )
    #fraction_distribution( input_path + '/data_fraction_distribution', N=10000, A=A, B=B, p_target=0.7, epsilon_p=0.03 )
    #fraction_distribution( input_path + '/data_fraction_distribution', N=100000, A=A, B=B, p_target=0.7, epsilon_p=0.03 )
    #fraction_distribution( input_path + '/data_fraction_distribution', N=1000000, A=A, B=B, p_target=0.7, epsilon_p=0.03 )


    A = 0.125
    B = 0.25
    #fraction_distribution( input_path + '/data_fraction_distribution', N=1000000, A=A, B=B, p_target=0.7, epsilon_p=0.05 )

    plt.show()

