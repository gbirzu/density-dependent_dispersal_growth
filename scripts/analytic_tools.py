import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.stats as stats
import glob


# Define variables
output_path = '../plots/'

def strip(s): #separates numbers from letters in string
    head = s.strip('-.0123456789')
    tail = s[len(head):]
    return head, tail

def get_variables(name):
    het_name = name.split('/')[-1].split('.')#Get file name
    if 'txt' in het_name:
        het_name.remove('txt') #remove ending
    het_name = '.'.join(het_name) #reform name
    aux = [strip(s) for s in het_name.split('_')]
    #default values if none found
    g = 1.0
    N = 1.0
    m0 = 1.0
    m1 = 1.0
    x = 1.0
    for s in aux:
        if s[0] == 'gf':
            g = float(s[1])
        elif s[0] == 'N':
            N = int(s[1])
        #elif s[0] == 'm':
        #    print s[1]
        #    [head, tail] = s[1].split('-')
        #    if head == '0':
        #        m0 = float(tail)
        #    elif head == '1':
        #        m1 = float(tail)
        elif s[0] == 'x':
            x = int(s[1])
    return g, N, m0, m1, x

def vA( r0, D0, A ):
    '''Computes exact velocity of wave with
    D(n) = D0*(1+An/N) and r(n) = r0*n*(1-n/N).
    Analytical result from Kawasaki et al., Theor. Ecol. (2017)
    '''
    if A >= 2.0:
        v = np.sqrt(r0 * D0) * (np.sqrt(A/2) + np.sqrt(2/A))
    else:
        v = 2 * np.sqrt(r0 * D0)
    return v


def vB( r0, D, B ):
    if B >= 2.0:
        v = np.sqrt(r0*D*B/2.)*(1 + 2./B)
    else:
        v = 2.*np.sqrt(r0*D)
    return v

def pqk_model(x, x_0, m_q, r_0, p, q):
    # Computes velocity and profile for exactly solvable models defined in SI
    D_q = m_q # factor of 2 from converting to hermitian form
    lambda_q = np.sqrt( ( q + 1.0 )*D_q/r_0 )
    #v = np.sqrt(dn*r0/(p + q))

    v = np.sqrt( ( q + 1.0 )/( p + q ) )*D_q/lambda_q

    prof = 1. - np.exp( ( x - x_0 )/lambda_q )
    prof[np.where(prof < 0.)[0]] = 0.
    prof = prof**(1./q)

    return v, prof

def theory_profile(z, z0, q):
    v = np.sqrt(0.25*0.01/(1 + q))
    return (1. - np.exp(v*(z - z0)/(0.25*np.sqrt(q + 1))))**(1./q)

def per_capita_growth( rho, r0, B ):
    return r0*( 1. - rho )*( 1. + B*rho )

def dispersal_rate( rho, D0, A, delta ):
    return D0*( 1. + A*rho**delta )

def velocity():
    vel_data = np.load('../data/linear_data/velocity_N1000000_gf0.01_m0-0.0_m1-0.25_avg_demes300.npy')
    l = len(vel_data)
    t_arr = vel_data.T[0]
    pop_arr = vel_data.T[1]

    slope, intercept, r_value, p_value, std_err = stats.linregress(t_arr[l/2:], pop_arr[l/2:])
    return slope

def ancestry_z(z, v, x_arr, prof_arr):
    i_z = np.argmin(abs(x_arr - z))
    if prof_arr[i_z] == 0.:
        ancestry = 0.
    else:
        x_int = x_arr[np.nonzero(prof_arr)]
        prof_int = prof_arr[np.nonzero(prof_arr)]
        prof_int = prof_arr[np.where(x_int <= z)]
        x_int = x_arr[np.where(x_int <= z)]

        integral = (v/0.25)*integrate.simps(1./prof_int, x_int)
        ancestry = prof_arr[i_z]**2*np.exp(integral)
    return ancestry

def ancestry_numeric(profile_array, v, m0, A, x_start=4):
    l_arr = len(profile_array)
    x_array = range( l_arr )
    r0 = 0.01
    D0 = m0/2.

    # Find index w/ last non-zero density
    i_max = np.where( profile_array == 0. )[ 0 ][ 0 ] - 1

    ancestry = np.zeros( l_arr )

    # Get first non-zero integral as baseline
    zetamax_index = i_max - x_start
    baseline_integral = ( v/D0 )*integrate.trapz( 1./( 1.0 + A*profile_array[ zetamax_index:i_max ] ), x_array[ zetamax_index:i_max ] )
    print baseline_integral
    for i in range( x_start, i_max ):
        zeta_index = i_max - i
        integration_profile = profile_array[ zeta_index:i_max ]
        integration_x = x_array[ zeta_index:i_max ]
        integral = ( v/D0 )*integrate.trapz( 1./( 1.0 + A*integration_profile ), integration_x )
        #if integral < 0.0:
            #integral -= baseline_integral
        #print 1./prof_integration[ :i ], x_integration[ :i ]
        ancestry[ zeta_index ] = ( profile_array[ zeta_index ]**2 )*np.exp( -integral )
    ancestry /= integrate.trapz( ancestry[ :zetamax_index ], x_array[ :zetamax_index ] )

    return np.array( [ x_array, ancestry ] )

def ancestry_analytical(x, x0, r0, m1):
    D1 = m1/4.0
    lm1 = np.sqrt(2*D1/r0) # Profile decay const
    ancestry = (2./lm1)*np.exp((x - x0)/lm1)*(1. - np.exp((x - x0)/lm1))
    ancestry[np.where(x > x0)[0]] = 0.
    return ancestry

def fixation_analytical(x, x0, r0, m1):
    D1 = m1/2.
    lm1 = np.sqrt(2*D1/r0) # Profile decay const
    fixation = np.exp((x - x0)/lm1)/lm1
    fixation[np.where(x > x0)[0]] = 0.
    return fixation

def exact_exponent( vr ):
    alpha = 2.*np.sqrt( 1. - 1./vr**2 )/( 1 - np.sqrt( 1. - 1./vr**2 ) )
    if ( type( vr ) is np.ndarray ) or ( isinstance( vr, list ) ):
        alpha[ np.where( alpha > 1.0 )[ 0 ] ] = 1.
    else:
        if alpha > 1:
            alpha = 1
    return alpha

def allee_profile( gf, migr, fstr, x ):
    D = migr/2.
    if fstr > -0.5:
        prof = 1./( 1. + np.exp( np.sqrt( gf/( 2.*D ) )*x ) )
    else:
        prof = 1./( 1. + np.exp( np.sqrt( gf*abs( fstr )/D )*x ) )
    return prof

def gaussian_pdf( x, x0, sigma ):
    return ( 1./np.sqrt( 2*np.pi*sigma**2 ) )*np.exp( -( x - x0 )**2/( 2.*sigma**2 ) )

def linear_cooperativity( A, B ):
    '''Takes A and B for linear dispersal model and returns:
        0 if wave is pulled
        1 if wave is semi-pushed
        2 if wave is fully-pushed'''
    if 2*A + B <= 2:
        return 0
    elif 2*A + B < 4:
        return 1
    else:
        return 2

def plot_ancestry():
    prof_file = '../data/linear_data/profile_N1000000_gf0.01_m0-0.0_m1-0.25_avg_demes300.npy'
    vel_file = '../data/linear_data/velocity_N1000000_gf0.01_m0-0.0_m1-0.25_avg_demes300.npy'
    r0, N, m0, m1, x = get_variables(prof_file)

    x_arr = np.arange(0, 300, 0.01)
    x_fix = np.arange(130, 180, 0.1)
    x0 = 157.
    p = 1.
    q = 1.
    v, prof_arr = pqk_model(x_arr, x0, 0.25, 0.01, q, p)
    v, prof_fix = pqk_model(x_fix, x0, 0.25, 0.01, q, p)

    prof_data = np.load(prof_file)
    vel_data = np.load(vel_file)
    x_sim = prof_data.T[0]
    prof_sim = (prof_data.T[1] + prof_data.T[2])/N
    #ancestry_num = ancestry_numeric(x_arr, prof_arr, v)
    ancestry_fix = ancestry_numeric(x_fix, prof_fix, v)
    #ancestry_num = ancestry_numeric(x_arr, prof_arr, v)
    i0 = np.argmin(abs(x_arr - x_fix[0]))
    i1 = np.argmin(abs(x_arr - x_fix[-1]))

    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_subplot(121)
    ax.set_xlabel('position [demes]')
    ax.set_xlim([130, 160])
    ax.set_ylim([0., 1.1])
    #ax.set_yscale('log')
    ax.plot(x_arr, prof_arr, c='b', label='density - theory')
    ax.scatter(x_sim, prof_sim, c='k', label='population density')
    #ax.plot(x_arr[i0:i1], ancestry_num[i0:i1], c='r', label='ancestry distribution - numerical')
    #ax.plot(x_fix, ancestry_fix, c='r', label='ancestry distribution - numerical')
    ax.plot(x_fix, 10*ancestry_analytical(x_fix, x0, r0, m1), ls='--', c='r', label='ancestry - simplified formula')
    ax.legend(loc='best', scatterpoints=1)

    ax = fig.add_subplot(122)
    ax.set_xlabel('time [generations]')
    ax.set_ylabel('distance travelled')
    ax.scatter(vel_data.T[0], vel_data.T[1], s=10, lw=0.5, facecolor='none', edgecolor='b', label='simulations')
    ax.plot(vel_data.T[0], v*vel_data.T[0] + 0.1*N, c='k', label='theory', lw=1.5)
    ax.legend(loc='best')

if __name__=='__main__':
    s = 3.0/( 2*np.sqrt( 2 ) ) - 1
    fractions = np.array( [ 1.0/4, 1.0/3, 1.0/2, 2.0/3, 3.0/4, 4.0/5 ] )
    nu_array = np.ones( len( fractions ) ) + s*fractions
    print fractions
    print nu_array
    print exact_exponent( nu_array )


    plt.show()

