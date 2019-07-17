#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

// Define RNG type
typedef gsl_rng Random;
typedef gsl_rng_type Random_type;

// Default parameters; used if none provided
const unsigned NUMBER_DEMES = 300;
const unsigned NUMBER_STRAINS = 2;

class stepping_stone{
    public:
        Random *rng;
        const unsigned demes;
        unsigned number_strains;
        std::vector< std::vector< int > > population;
        std::vector< std::vector< double > > density;

        //stepping_stone( int, unsigned, unsigned, double, bool );
        stepping_stone( int time=0,
                unsigned demes=NUMBER_DEMES,
                unsigned number_strains=NUMBER_STRAINS,
                double centering=0.4,
                bool stochastic=true );

    private:
        void initialize_rng();
        void initialize_arrays();
};
