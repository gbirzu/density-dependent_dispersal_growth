#include "stepping_stone.hpp"

/*stepping_stone::stepping_stone( int time=0,
        unsigned demes=NUMBER_DEMES,
        unsigned number_strains=NUMBER_STRAINS,
        double centering=0.4,
        bool stochastic=true ): demes( demes ), number_strains( number_strains ){
    initialize_arrays();
    if( stochastic == true ){
        initialize_rng();
    }
}*/

stepping_stone::stepping_stone( int time,
        unsigned demes,
        unsigned number_strains,
        double centering,
        bool stochastic ): demes( demes ), number_strains( number_strains ){
    initialize_arrays();
    if( stochastic == true ){
        initialize_rng();
    }
}


void stepping_stone::initialize_rng(){
    const Random_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    rng = gsl_rng_alloc(T);

    // Get Seed for PRG
    int sysRandom; // Seed variable
    std::ifstream firand;
    firand.open("/dev/urandom", std::ios_base::in | std::ios_base::binary); // Open urandom as binary file
    firand.read((char*) (&sysRandom), sizeof sysRandom); // Read from urandom
    firand.close();
    std::cout << "Random seed: " << sysRandom << std::endl;
    gsl_rng_set(rng, sysRandom); // Seed RNG
}

void stepping_stone::initialize_arrays(){
        population.resize( demes, std::vector< int >( number_strains ) );
        density.resize( demes, std::vector< double >( number_strains ) );
}

/*
// Define RNG type
typedef gsl_rng Random;
typedef gsl_rng_type Random_type;

// Default parameters; used if none provided
const unsigned NUMBER_DEMES = 300;
const unsigned NUMBER_STRAINS = 2;

class stepping_stone{
    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//
    // Class sets up 1D stepping-stone model for simulations of traveling fronts
     * Initializes RNG and sets up arrays for: population count, density, intermediate values
     * Assumes field is discreet;
     * population[][] stores the number of each strain; first index stores position, second strain number
     * density[][] is population[][]/carrying_capacity
    //
    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//
    
    public:
        Random *rng;
        const unsigned demes;
        unsigned number_strains;
        std::vector< std::vector< int > > population;
        std::vector< std::vector< double > > density;

        stepping_stone( int time=0,
                unsigned demes=NUMBER_DEMES,
                unsigned number_strains=NUMBER_STRAINS,
                double centering=0.4,
                bool stochastic=true ): demes( demes ), number_strains( number_strains ){
            initialize_arrays();
            if( stochastic == true ){
                initialize_rng();
            }
        }

    private:
        void initialize_rng(){
            const Random_type *T;
            gsl_rng_env_setup();
            T = gsl_rng_mt19937;
            rng = gsl_rng_alloc(T);

            // Get Seed for PRG
            int sysRandom; // Seed variable
            std::ifstream firand;
            firand.open("/dev/urandom", std::ios_base::in | std::ios_base::binary); // Open urandom as binary file
            firand.read((char*) (&sysRandom), sizeof sysRandom); // Read from urandom
            firand.close();
            std::cout << "Random seed: " << sysRandom << std::endl;
            gsl_rng_set(rng, sysRandom); // Seed RNG
        }

        void initialize_arrays(){
            population.resize( demes, std::vector< int >( number_strains ) );
            density.resize( demes, std::vector< double >( number_strains ) );
        }
};
*/
