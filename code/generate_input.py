import numpy as np
import glob

def change_extension( file_name, new_ext ):
    components = file_name.split( '.' )
    components[ -1 ] = new_ext
    new_name = '.'.join( components )
    return new_name


def change_head( file_name, new_head ):
    components = file_name.split( '_' )
    components[ 0 ] = new_head
    new_name = '_'.join( components )
    return new_name


def convert_to_txt( input_file, output_file ):
    # Takes .npy profile and converts it to .txt
    # only outputs total density
    profile_input = np.load( input_file )
    input_array = ( profile_input.T[ :2 ] ).T

    # Add strain numbers
    input_array.T[ 1 ] = input_array.T[ 1 ] + profile_input.T[ 2 ]

    input_array.T[ 1 ] /= max( input_array.T[ 1 ] )
    input_array.T[ 0 ] = ( input_array.T[ 0 ] ).astype( int )
    np.savetxt( output_file, input_array, delimiter=',', fmt='%d,%1.4e' )

    # Save last non-zero element
    x_max = input_array[ np.where( input_array.T[ 1 ] == 0.0 )[ 0 ][ 0 ], 0 ] - 1
    f_xmax = open( change_head( output_file, "xmax" ), "w" )
    f_xmax.write( str( int( x_max ) ) )
    f_xmax.close()


def convert_all_files():
    # Searches for all profile files and converts them to .txt
    profile_npy_list = glob.glob( "profile_*.npy" )
    for profile_npy in profile_npy_list:
        output_profile = change_extension( profile_npy, 'txt' )
        convert_to_txt( profile_npy, output_profile )

if __name__ == '__main__':
    convert_all_files()
