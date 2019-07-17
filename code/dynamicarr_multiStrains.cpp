#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <memory.h> //needed for memcpy
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

//define default parameters of simulation
const int Lx = 400; // Number of demes in x direction; must be even
const int Ly = 400; // Number of demes in y direction; must be even
const int NStrains = 10; // Number of strains
const int NDeme = Lx*Ly; // Number of demes in box
const int ArraySize = NDeme*NStrains; // Size of density array
unsigned int NGen = 5; // Number of generations to run
unsigned long K = 100; // Size of deme
double Gf = 0.001; // Growth rate
double Fstar = -1.0; // f*
double M = 0.08; // Migration rate
int R0 = 8; // Radius of colony at t=0

unsigned nHetPoints = 100;
unsigned nProfPoints = 1;
int Run = 0;

int *deme; // Data array
int *newDeme; // Previous generation array

double calcHet() //deme average
{
    double h = 0.0;
    int counter = 0;
    for (int i = 0; i < NDeme; i++)
    {
        double population = 0.0;
        for (int j = 0; j < NStrains; j++)
            population += double(deme[i + j*NDeme]);
        if (population > 0.0)
        {
            for (int j = 0; j < NStrains; j++)
                h += (deme[i + j*NDeme]/population)*(1.0 - deme[i + j*NDeme]/population);
            counter++;
        }
    }
    return h/counter;
}

gsl_rng *initializeRng()
{
    gsl_rng *rng;
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    rng = gsl_rng_alloc(T);

    // Get Seed for PRG
    int sysRandom; // Seed variable
    ifstream firand;
    firand.open("/dev/urandom", ios_base::in | ios_base::binary); // Open urandom as binary file
    firand.read((char*) (&sysRandom), sizeof sysRandom); // Read from urandom
    firand.close();
    cout << "Random seed: " << sysRandom << endl;
    gsl_rng_set(rng, sysRandom); // Seed RNG
    return rng;
}


void initializeArrays(gsl_rng *r)
{
    deme = new int [ArraySize]; 
    newDeme = new int [ArraySize];

    // Initialize state
    double initialFraction = 1./double(NStrains);
    int xCenter = Lx/2;
    int yCenter = Ly/2;
    for (int x = 0; x < Lx; x++)
        for (int y = 0; y < Ly; y++)
        {
            int i = x + y*Lx; // Index of (x,y)
            int dx = x - xCenter; // Delta x from center of box
            int dy = y - yCenter; // Delta y from center of box
            double pAux[NStrains + 1]; // Strain probabilities
            unsigned int nAux[NStrains + 1]; // Strain numbers

            if (pow(dx, 2) + pow(dy, 2) <= pow(R0, 2))
            {
                // If (x,y) is inside radius of the initial colony, assign strains from uniform distribution
                for (int j = 0; j < NStrains; j++)
                    pAux[j] = initialFraction; // species probabilities
                pAux[NStrains] = 0.0;
                gsl_ran_multinomial(r, NStrains + 1, K, pAux, nAux); // get next generation

                for (int j = 0; j < NStrains; j++)
                    deme[i + j*NDeme] = nAux[j];
            }
            else
                // Else deme is empty
                for (int j = 0; j < NStrains; j++)
                {
                    deme[i + j*NDeme] = 0;
                    newDeme[i + j*NDeme] = 0; // Always initialize newDeme to zero
                }
        }

}

void initializeOutputVariables(int *tHetInterval, int *tProfInterval)
{
    if (NGen > nHetPoints)
        *tHetInterval = int(NGen/nHetPoints); // Number of generation per data save
    else
        *tHetInterval = 1; // Number of generation per data save

    if (NGen > nProfPoints)
        *tProfInterval = int(NGen/nProfPoints); // Number of generation per data save
    else
        *tProfInterval = 1; // Number of generation per data save
}

void openOutputFiles(ofstream& fHet, ofstream& fStrains)
{
    // Convert paramenters to strings
    ostringstream strK, strFstar, strG, strM, strRun, strNDeme, strNStrains, strNGen, strR0;
    strK << K;
    strFstar << setprecision(4) << Fstar;
    strG << Gf;
    strM << M;
    strRun << Run;
    strNDeme << NDeme;
    strNStrains << NStrains;
    strNGen << NGen;
    strR0 << R0;

    // Create file names
    string termination = "_t" + strNGen.str() + "_rinit" + strR0.str() + "_strains" + strNStrains.str() + ".txt";
    string hetName = "hetero_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_fstr" + strFstar.str() + "_run" + strRun.str() + termination;
    string strainsName = "strains_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_fstr" + strFstar.str() + "_run" + strRun.str() + termination;

    // Open files
    fHet.open(hetName);
    fStrains.open(strainsName);

    if (fHet.is_open())
        cout << "Heterozygosity ouput file succesfully opened!" << endl;
    else
        cout << "ERROR: Heterozygosity ouput file failed to open!" << endl;

    if (fStrains.is_open())
    {
        cout << "Profile ouput file succesfully opened!" << endl;
        //Print parameters to file
        fStrains << Lx << ',' << Lx << endl;
        fStrains << nProfPoints << endl;
    }
    else
        cout << "ERROR: Profile ouput file failed to open!" << endl;
}


double migrateFractions(double f, double fUp, double fDown, double fRight, double fLeft)
{
    double newFractions = (1 - M)*f + (M/4.0)*fUp + (M/4.0)*fDown + (M/4.0)*fRight + (M/4.0)*fLeft;//migrate
    return newFractions;
}

void growthFractions(double *newFractions)
{
    // Takes pointer to array of species fractions and updates them according to the growth function
    double wS = 1.0; // Species fitness
    double speciesFraction = 0.0;
    long double wAvg = 0.0; // Average fitness
    long double wV; // Vacancy fitness

    for (int i = 0; i < NStrains; i++)
        speciesFraction += newFractions[i];

    wV = 1.0 + Gf*Fstar - Gf*speciesFraction; // Assign wV
    wAvg = (1. - speciesFraction)*wV + speciesFraction*wS; // Assign wAvg
    // Calculate updated fractions
    for(int i = 0; i < NStrains; i++)
        newFractions[i] *= wS/wAvg;
}

void assignNextGeneration(int i, double *newFractions, gsl_rng *r)
{
    // Updates newDeme at i = x + y*Lx stochastically

    double pAux[NStrains + 1];
    unsigned int nAux[NStrains + 1];
    double speciesFraction = 0.0;

    /////////////// Stochastic assignment of species ///////////
    for (int j = 0; j < NStrains; j++)
    {
        pAux[j] = newFractions[j]; // Species probabilities
        speciesFraction += newFractions[j];
    }
    pAux[NStrains] = 1.0 - speciesFraction;// Vacancy probablity
    gsl_ran_multinomial(r, NStrains + 1, K, pAux, nAux); // Get next generation
    for (int j = 0; j < NStrains; j++)
        newDeme[i + j*NDeme] = nAux[j]; // Assign next generation
}

void updateBulk(gsl_rng *r)
{
    int i, i_up, i_down, i_right, i_left;
    double *newFractions;
    for (int x = 1; x < Lx - 1; x++)
        for (int y = 1; y < Ly - 1; y++)
        {
            newFractions = new double [NStrains];
            i = x + y*Lx;
            i_up = x + (y + 1)*Lx;
            i_down = x + (y - 1)*Lx;
            i_right = x + 1 + y*Lx;
            i_left = x - 1 + y*Lx;
            
            ///////////////// Calculate species fractions //////////////////
            //-------------- Migration -------------------//
            for (int j = 0; j < NStrains; j++)
                newFractions[j] = migrateFractions(deme[i + j*NDeme]/double(K), deme[i_up + j*NDeme]/double(K), deme[i_down + j*NDeme]/double(K), deme[i_left + j*NDeme]/double(K), deme[i_right + j*NDeme]/double(K));
            //-------------------------------------------------------------------------//

            //--------------------------- Growth -------------------//
            growthFractions(newFractions);
            //-------------------------------------------------------------------------//

            //------------------ Sample next generation stochastically -----------------//
            assignNextGeneration(i, newFractions, r);
            //-------------------------------------------------------------------------//

            delete newFractions;
        }
}

void applyBC()
{
    int i_b, i_bint, i_t, i_tint, i_l, i_lint, i_r, i_rint;
    // Updates edges of newDeme using reflecting BC
    for (int x = 0; x < Lx; x++)
    {
        i_b = x; // Index of bottom boundary
        i_bint = x + Lx; // Index of interior of the bottom boundary
        i_t = x + (Ly - 1)*Lx; // Index of bottom boundary
        i_tint = x + (Ly - 2)*Lx; // Index of interior of the bottom boundary

        // Apply reflecting BC
        for (int j = 0; j < NStrains; j++)
        {
            newDeme[i_b] = newDeme[i_bint];
            newDeme[i_t] = newDeme[i_tint];
        }
    }

    for (int y = 0; y < Ly; y++)
    {
        i_l = y*Lx;
        i_lint = 1 + y*Lx;
        i_r = Lx - 1 + y*Lx;
        i_rint = Lx - 2 + y*Lx;

        // Apply reflecting BC
        for (int j = 0; j < NStrains; j++)
        {
            newDeme[i_l] = newDeme[i_lint];
            newDeme[i_r] = newDeme[i_rint];
        }
    }
}

void saveProfile(int t, ofstream& fStrains, bool end)
{
    fStrains << t << endl;
    for (int i = 0; i < NDeme; i++)
    {
        fStrains << i;
        for (int j = 0; j < NStrains; j++)
            fStrains << ',' << deme[i + j*NDeme]/double(K);
        fStrains << endl;
    }
    if (end == false)
        fStrains << endl;
}


int main(int argc, char *argv[])
{
    clock_t c_i = clock(); // Start program timer
    int c;

    // Read parameters from input flags
    while ((c = getopt (argc, argv, "r:f:m:g:K:t:i:")) != -1)
    {
        if (c == 'r')
            Run  = atoi(optarg);
        else if (c == 'f')
            Fstar = atof(optarg);
        else if (c == 'm')
            M = atof(optarg);
        else if (c == 'g')
            Gf = atof(optarg);
        else if (c == 'K')
            K = atoi(optarg);
        else if (c == 't')
            NGen = atoi(optarg);
        else if (c == 'i')
            R0 = atof(optarg);
    }

    cout << "Run variable: " << Run << '\t' << "Fstar var: " << Fstar << "\t Migr var: " << M << "\t gf: " << Gf << "\t K: " << K << "\t Generations: " << NGen << "\t Number of strains: " << NStrains << endl;

    // Declare simulation variables
    unsigned t = 0; // Time variable
    int tHetInterval, tProfInterval;
    ofstream fHet, fStrains;
    vector<double> pop; // Population vector
    vector<double> het; // Heterozygosity vector

    // Initialization 
    gsl_rng *rng = initializeRng();
    initializeArrays(rng);
    openOutputFiles(fHet, fStrains);
    initializeOutputVariables(&tHetInterval, &tProfInterval);

    cout << endl << "Initial heterozygosity: " << calcHet() << endl;

    //------------------------------------------------------------------------//
    for (; t < NGen; t++)
    {
        copy(deme, deme + ArraySize, newDeme);
        updateBulk(rng);
        applyBC();
        copy(newDeme, newDeme + ArraySize, deme);

        // ------------------------- Store data -------------------------//
        if (t%tHetInterval == 0)
            het.push_back(calcHet());//save heterozygosity

        //Profile file
        if (t%tProfInterval == 0)
            saveProfile(t, fStrains, true); // Save final profile
        //------------------------------------------------------------------------//
    }
    delete deme;
    delete newDeme;
    clock_t c_f = clock();

    // ------------------------- Write output to files -------------------------//
    //------------------------------------------------------------------------//
    saveProfile(t, fStrains, true); // Save final profile

    //Heterozygosity file
    for (unsigned int i = 0; i < het.size(); i++)
       fHet << tHetInterval*i << ',' << het[i] << endl;
    //------------------------------------------------------------------------//
   
    cout << "Run time: " << double(c_f - c_i)/CLOCKS_PER_SEC << endl;
    cout << endl << "Final heterozygosity: " << calcHet() << endl;
    return 0;
}
