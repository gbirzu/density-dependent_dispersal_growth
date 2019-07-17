#include "stepping_stone.hpp"

int main( int argc, char *argv[] ){
    stepping_stone lattice( 0, NUMBER_DEMES, NUMBER_STRAINS, 0.5, true );
    std::cout << lattice.density[ 0 ][ 0 ] << ' ' << lattice.density[ 0 ][ 1 ] << std::endl;
    std::cout << "New line" << std::endl;
    return 0;
}
