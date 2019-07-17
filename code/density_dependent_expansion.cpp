/*-----------------------------------------------------------------------*/
/* This file contains a class that sets up 1D stepping-stone model for
 * simulations of traveling fronts
*/
/*-----------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

// For random number generators
typedef gsl_rng Random;
typedef gsl_rng_type Random_type;

// Define constants
const unsigned NDEMES = 300;

// Class definitions
class environment
{
    // Initializes RNG and sets up arrays for: population count, density, intermediate values
    public:
        const unsigned demes;
        unsigned array_size;
        int* density;
        double* rho;
        Random *rng;
        std::ofstream f_het, f_profile, f_total_density;

        environment(int time=0, unsigned demes=NDEMES, double centering=0.4): demes(demes)
        {
            array_size = 2*demes;
            initialize_rng();
            initialize_arrays();
        }
    private:

        void initialize_rng()
        {
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

        void initialize_arrays() 
        {
            density = new int [array_size]; 
            rho = new double [array_size];

            for (int i = 0; i < array_size; i++)
            {
                density[i] = 0.0;
                rho[i] = 0.0;
            }
        }
};

class population
{
    // Runs range expansion simulations and tracks velocity, profile, heterozygosity, etc.
    const char* migration_type;
    public:
        environment env;
        double migration;
        double migration_cooperativity;
        double growth;
        double growth_cooperativity;
        const unsigned carrying_capacity;
        const int simulation_time;
        int demes_shifted;

        std::vector<double> het_vector;
        std::vector<double> distance_vector;
        std::vector<std::vector<double>> final_rho;

        population(const char* migration_type="linear", double migration=0.25, double migration_cooperativity=1.0, double growth=0.1, double growth_cooperativity=1.0, unsigned carrying_capacity=1000, const int simulation_time=0, bool equilibrate=true): migration_type(migration_type), migration(migration), migration_cooperativity(migration_cooperativity), growth(growth), growth_cooperativity(growth_cooperativity), carrying_capacity(carrying_capacity), simulation_time(simulation_time)
        {
            demes_shifted = 0;
            if(simulation_time > 0)
            {
                het_vector.reserve(simulation_time);
                distance_vector.reserve(simulation_time);
            }
            else
            {
                het_vector.reserve(10000000);
                distance_vector.reserve(10000000);
            }
        }

        double get_heterozygosity()
        {
            double h = 0.0;
            int counter = 0;

            for (int x = 0; x < env.demes; x++)
            {
                double population = double(env.density[x] + env.density[x + env.demes]);
                if (population > 0.0)
                {
                    h += 2.*(env.density[x]/population)*(env.density[x + env.demes]/population);
                    counter++;
                }
            }
            return h/counter;
        }


        double get_distance()
        {
            double distance = 0.0;
            for (int x = 0; x < env.demes; x++)
                distance += env.rho[x] + env.rho[x + env.demes];
            return distance;
        }


        void simulation_step()
        {
            update_fractions();
            migrate();
            grow();
            assign_densities();

            // Update data
            het_vector.push_back(get_heterozygosity());
            distance_vector.push_back(get_distance());
        }

        void simulate()
        {
            // Runs simulation for simulation_time and saves results;
            // if simulation_time= 0 simulations runs until heterozygosity = 0
            assign_initial_conditions();
            het_vector.push_back(get_heterozygosity());
            distance_vector.push_back(get_distance());
            if(simulation_time == 0)
            {
                int time = 0;
                while(het_vector.back() > 0.0)
                {
                    simulation_step();
                    time++;
                }

                // Save final profile
                for(int x = 0; x < env.demes; x++)
                {
                    std::vector<double> rho_x(env.rho[x], env.rho[x + env.demes]);
                    final_rho.push_back(rho_x);
                }

            }
        }


        void output_data(int run_number)
        {
        }

    private:
        void assign_initial_conditions()
        {
            // Initialize state
            double initialFraction = 0.5;
            for (int x = 0; x < env.demes; x++)
            {
                double pAux[3]; // Strain probabilities
                unsigned int nAux[3]; // Strain numbers

                // Fill first half of simulation box
                if (x < env.demes/2)
                {
                    pAux[0] = initialFraction; // species probabilities
                    pAux[1] = initialFraction; // species probabilities
                    pAux[2] = 0.0;
                    gsl_ran_multinomial(env.rng, 3, carrying_capacity, pAux, nAux); // get next generation
                    env.density[x] = nAux[0];
                    env.density[x + env.demes] = nAux[1];
                }
                else
                {
                    // Else deme is empty
                    env.density[x] = 0;
                    env.density[x + env.demes] = 0;
                }
                env.rho[x] = env.density[x]/double(carrying_capacity);
                env.rho[x + env.demes] = env.density[x + env.demes]/double(carrying_capacity);
            }
        }
        

        double mr_linear(double rho)
        {
            return migration*(1. + migration_cooperativity*rho);
        }

        double mr_sqroot(double rho)
        {
            return migration*(1. + migration_cooperativity*sqrt(rho));
        }

        double mr_quadratic(double rho)
        {
            return migration*(1. + migration_cooperativity*pow(rho, 2));
        }


        std::vector<double> get_migration_rates()
        {
            std::vector<double> mr_rate(env.demes);
            std::vector<double> rho_total(env.demes);
            for(int i = 0; i < env.demes; i++)
            {
                rho_total[i] = (env.rho[i] + env.rho[i + env.demes]);
                if(strcmp(migration_type, "linear") == 0)
                    mr_rate[i] = mr_linear(rho_total[i]);
                else if (strcmp(migration_type, "quadratic") == 0)
                    mr_rate[i] = mr_quadratic(rho_total[i]);
                else if (strcmp(migration_type, "square_root") == 0)
                    mr_rate[i] = mr_sqroot(rho_total[i]);
            }
            return mr_rate;
        }


        void update_fractions()
        {
            for(int i = 0; i < env.demes; i++)
            {
                env.rho[i] = env.density[i]/double(carrying_capacity);
                env.rho[i + env.demes] = env.density[i + env.demes]/double(carrying_capacity);
            }
        }


        void migrate()
        {
            // Migrates rhos w/ BC
            double* new_rho = new double [env.array_size];
            std::vector<double> mr_rate = get_migration_rates();
            for(int x = 1; x < env.demes - 1; x++)
            {
                new_rho[x] = (1. - mr_rate[x])*env.rho[x] + 0.5*mr_rate[x - 1]*env.rho[x - 1] + 0.5*mr_rate[x + 1]*env.rho[x + 1];
                new_rho[x + env.demes] = (1. - mr_rate[x])*env.rho[x + env.demes] + 0.5*mr_rate[x - 1]*env.rho[x - 1 + env.demes] + 0.5*mr_rate[x + 1]*env.rho[x + 1 + env.demes];
            }

            // Apply reflective boundary conditions
            // Left boundary
            new_rho[0] = (1. - 0.5*mr_rate[0])*env.rho[0] + 0.5*mr_rate[1]*env.rho[1];
            new_rho[env.demes] = (1. - 0.5*mr_rate[0])*env.rho[env.demes] + 0.5*mr_rate[1]*env.rho[env.demes + 1];
 
            // Right boundary
            new_rho[env.demes - 1] = (1. - 0.5*mr_rate[env.demes - 1])*env.rho[env.demes - 1] + 0.5*mr_rate[env.demes - 2]*env.rho[env.demes - 2];
            new_rho[2*env.demes - 1] = (1. - 0.5*mr_rate[env.demes - 1])*env.rho[2*env.demes - 1] + 0.5*mr_rate[env.demes - 2]*env.rho[2*env.demes - 2];
 
            std::copy(new_rho, new_rho + env.array_size, env.rho);
        }


        void grow()
        {
            // Add growth to rho
            double rho_total = 0.;
            double growth_rate = 0.;
            for(int x = 0; x < env.demes; x++)
            {
                rho_total = env.rho[x] + env.rho[x + env.demes]; 
                growth_rate = growth*(1. - rho_total)*(1. + growth_cooperativity*rho_total);
                env.rho[x] += growth_rate*env.rho[x];
                env.rho[x + env.demes] += growth_rate*env.rho[x + env.demes];
                if (env.rho[x] < 0.)
                    env.rho[x] = 0.0;
                if (env.rho[x + env.demes] < 0.)
                    env.rho[x + env.demes] = 0.0;
            }
        }


        void assign_densities()
        {
            // Draw next generation from multinomial w/ probability rho
            double pAux[3];
            unsigned int nAux[3];

            for(int i = 0; i < env.demes; i++)
            {
                /////////////// Stochastic assignment of species ///////////
                pAux[0] = env.rho[i]; // Species probabilities
                pAux[1] = env.rho[i + env.demes]; // Species probabilities
                if (pAux[0] + pAux[1] == 1.)
                    pAux[2] = 0.;
                else
                    pAux[2] = 1.0 - env.rho[i] - env.rho[i + env.demes];// Vacancy probablity
                gsl_ran_multinomial(env.rng, 3, carrying_capacity, pAux, nAux); // Get next generation
                env.density[i] = nAux[0]; // Assign next generation
                env.density[i + env.demes] = nAux[1]; // Assign next generation
            }
        }


        std::string file_name(std::string head)
        {
            // Create file name with all relevant parameters; name begins w/ head
            std::ostringstream fname, termination;
            termination << "_t" << simulation_time << ".txt"; // Optional variables and file extension
            fname << head << "_N" << carrying_capacity << "_r0" << growth << "_m0" << migration << "_B" << growth_cooperativity << termination.str();
            return fname.str();
        }
};


// ################################################################## // 


// ------------------------------------------------------------------- //
int main(int argc, char* argv[]) 
{
    // Main function; reads in parameters, runs simulation and outputs heterozygosity, total population and final profile to .txt files
    if(argc < 9)
    {
        // Tell user how to use program
        std::cout << "Parameters input: " << argv[0] << " migration_type m0 A r0 B carrying_capacity time run_number " << std::endl;
        return 1;
    }
    else
    {
        // Read in parameters
        const char* migration_type = argv[1];
        double m0 = atof(argv[2]);
        double A = atof(argv[3]);
        double r0 = atof(argv[4]);
        double B = atof(argv[5]);
        unsigned N = atoi(argv[6]);
        unsigned t_end = atoi(argv[7]);
        unsigned run_number = atoi(argv[8]);

        population n(migration_type, m0, A, r0, B, N, t_end);
        n.simulate();
        n.output_data(run_number);
    }
    return 0;
}

// ------------------------------------------------------------------- //
