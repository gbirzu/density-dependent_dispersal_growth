// Code Version 6.1 C++
// Simulates expansions with linear dispersal and cubic growth
//  In addition to standard outputfiles, outputs the full profile at
//  100 evenly spaced timepoints

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

using namespace std ;

// Variables
const double pi = 3.14159 ;
const int nDeme = 300 ;
const int nSp = 2 ;
unsigned int nGen = 1000 ; 

// Growth and migration parameters
double r0 = 0.01 ; // Zero density growth rate
double B = 0.0 ; // Growth cooperativity
double m0 = 0.1 ; // Zero density migration rate
double A = 0.0 ; // Migration cooperativity
double delta ;     // Non-linear migration exponent
int populationSize = 1000 ;  // Carrying capacity
int xFixation = int(nDeme/4) ; // Default; check fixation at quarter of box
int batchNumber = 0 ;
int pDisp ;
int runNumber ;          // Run Number 

int nHetPoints = 1000 ; 
int nProfPoints = 100 ; 
int tHetInterval = int( nGen / nHetPoints) ; 
int tProfInterval = int( nGen / nProfPoints );
int widthCounter = 0 ;  // 
// ~~~~~~~~~~~~~~~~ //

// Flag Variables
int instDeme , maxDeme ; 
unsigned int deterministicUpdate = 0 ; // 0 for deterministic update; 1 for stochastic update
unsigned int checkLastDeme = 1; 
unsigned int widthCheckFlag = 0 ; 
// ~~~~~~~~~~~~~~~~~~~~~ //


// Array
long double sA[nDeme][nSp], sAP[nDeme][nSp], sAPTotal[nDeme] ;
long int dA[nDeme][nSp] ;

vector<double> pop ; 
vector<double> hetGlobal ; 
vector<double> fGlobal ; 
// ~~~~~~~~~~~~~ //


// RND GSL
const gsl_rng_type * T ; 
gsl_rng * r ; 
// ~~~~~~~~~~~~~~ //



// ------------------------------------------------------------------- //
// FUNCTIONS
// ------------------------------------------------------------------- //
void Initialize()
{
  for (int i = 0; i < nDeme ; i++ ) 
    {
      dA[i][0] = 0 ; 
      dA[i][1] = 0 ; 

      sA[i][0] = 0.0 ; 
      sA[i][1] = 0.0 ; 
    }

  for (int i = 0; i < nDeme / 2 ; i ++)
    {
      dA[i][0] = gsl_ran_binomial(r, 0.5, populationSize) ;
      dA[i][1] = populationSize - dA[i][0] ;

      sA[i][0] = dA[i][0] / float(populationSize) ;
      sA[i][1] = dA[i][1] / float(populationSize) ;
    }

  pDisp = 0 ; 
}
// ------------------------------------------------------------------- //


int fixationInitialize( int x, string finput_name )
{
    // Initialize all values
    for (int i = 0; i < nDeme ; i++ ) 
    {
        dA[i][0] = 0 ; 
        dA[i][1] = 0 ; 

        sA[i][0] = 0.0 ; 
        sA[i][1] = 0.0 ; 
    }

    // Read profile from input file
    ifstream finput( finput_name.c_str() ); 

    if( !finput )
      {
        cerr << "ERROR! Unable to open input file: " << finput_name << endl ;
        return -1 ;
      }

    string value_x, value_rho ;
    while ( !finput.eof() )
        {
        getline( finput, value_x, ',' ) ;
        getline( finput, value_rho, '\n' ) ;

        if (value_x != "")
            {
            if ( stoi(value_x) >= 0 and stoi(value_x) < nDeme )
                {
                if ( stoi(value_x) != x )
                    {
                    dA[stoi(value_x)][0] = long(populationSize*stod(value_rho)) ;
                    dA[stoi(value_x)][1] = 0 ;

                    sA[stoi(value_x)][0] = stod(value_rho) ;
                    sA[stoi(value_x)][1] = 0.0 ;
                    }
                else
                    {
                    dA[stoi(value_x)][0] = 0 ;
                    dA[stoi(value_x)][1] = long(populationSize*stod(value_rho)) ;

                    sA[stoi(value_x)][0] = 0.0 ;
                    sA[stoi(value_x)][1] = stod(value_rho) ;
                    }
                }
            }
        } 
    finput.close() ; 
    pDisp = 0 ; 
    return 0 ;
}

// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
void calcNGen( double pushedBoundary, double fullyPushedBoundary )
{
  // Uses heuristic to approximate simulation time
  // The boundaries between the phases are set by the inputs
  // Works for linear migration model only!
  
  double tBSC ;
  tBSC = (pow(log(populationSize),3))/(2*r0*pow(pi, 2)) ;
  if ( 2*A + B < pushedBoundary )
  {
      // Approximate phase boundary by straight line 
      nGen = max( int( 20*tBSC ), 50000 ) ;
  }
  else if ( 2*A + B >= pushedBoundary and 2*A + B < fullyPushedBoundary )
  {
      // Heuristic guess for simulation time; not the best, but kind of works 
      double exponent = ( 2*A + B  - pushedBoundary )/( fullyPushedBoundary - pushedBoundary ) ; // Linear approximation of alpha
      double tSemipushed = sqrt( m0/( 2*r0 ) )*pow( populationSize, exponent ) ; // Estimated semi-pushed simulation time

      nGen = max( 2*populationSize, int( 20*tSemipushed ) ); 
      // Ensure simulation time is at least 10000 generations
      nGen = max( int( nGen ), 50000 ) ; 
  }
  else
  {
      nGen = max( 20*populationSize, 50000 ) ;
  }

  tHetInterval = int( nGen / nHetPoints) ; 
  tProfInterval = int( nGen / nProfPoints ) ; 
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
double sumArr( long  arr[][nSp], int populationSize) // Only for dA
{
  double sum = 0.0 ; 
  for (int i = 0 ; i < nDeme ; i++)
    {
      sum += arr[i][0] / double(populationSize) + arr[i][1] / double(populationSize) ; 
    }
  return sum; 
}
// ------------------------------------------------------------------ //


// ------------------------------------------------------------------- //
double mig_rate( double sp_frac) // defines function for migration rate; uses m0 and A as parameters
{
    return m0*( 1.0 + A*sp_frac ) ;
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
void Mig_del_1() // D_1 term is polynomial of order 1. For Faster simulations 
{

  for (int i = 0; i < nDeme; i++) 
    {
      // Converts into Fractions
      sA[i][0] = double(dA[i][0]) / double(populationSize) ;
      sA[i][1] = double(dA[i][1]) / double(populationSize) ;

      // Specie Array Previous 
      sAP[i][0] = sA[i][0] ;
      sAP[i][1] = sA[i][1] ;      
      sAPTotal[i] = sAP[i][0] + sAP[i][1] ;
    }

  for (int i = 0; i < nDeme ; i ++ )
  {
      if ( i == 0 )
      {
          sA[i][0] = sAP[i][0] - 0.5 * mig_rate(sAPTotal[i]) * sAP[i][0] 
	        + 0.5 * mig_rate(sAPTotal[i+1]) * sAP[i+1][0] ;

          sA[i][1] = sAP[i][1] - 0.5 * mig_rate(sAPTotal[i]) * sAP[i][1] 
	        + 0.5 * mig_rate(sAPTotal[i+1]) * sAP[i+1][1] ;
      }

      else if ( i == nDeme - 1)
      {
          sA[i][0] = sAP[i][0] - 0.5 * mig_rate(sAPTotal[i]) * sAP[i][0] 
	        + 0.5 * mig_rate(sAPTotal[i-1]) * sAP[i-1][0] ;

          sA[i][1] = sAP[i][1] - 0.5 * mig_rate(sAPTotal[i]) * sAP[i][1] 
	        + 0.5 * mig_rate(sAPTotal[i-1]) * sAP[i-1][1] ;
      }

      else if ( i > 0 && i < nDeme ) 
      {
          sA[i][0] = sAP[i][0] - mig_rate(sAPTotal[i]) * sAP[i][0] 
              + 0.5 * mig_rate(sAPTotal[i-1]) * sAP[i-1][0]
              + 0.5 * mig_rate(sAPTotal[i+1]) * sAP[i+1][0] ;

          sA[i][1] = sAP[i][1] - mig_rate(sAPTotal[i]) * sAP[i][1] 
              + 0.5 * mig_rate(sAPTotal[i-1]) * sAP[i-1][1]
              + 0.5 * mig_rate(sAPTotal[i+1]) * sAP[i+1][1] ;
	  }
  }
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
void Growth()
{
  for (int i = 0; i<nDeme; i++) 
    {
      // Assign Prev. Array
      sAP[i][0] = sA[i][0] ;
      sAP[i][1] = sA[i][1] ;
      sAPTotal[i] = sAP[i][0] + sAP[i][1] ;
      

      // Implement Growth; per capita growth depends only on total number since species are neutral
      sA[i][0] = sAP[i][0] + r0 * sAP[i][0] * ( 1.0 - sAPTotal[i] ) * ( 1.0 + B * sAPTotal[i] );
      sA[i][1] = sAP[i][1] + r0 * sAP[i][1] * ( 1.0 - sAPTotal[i] ) * (1.0 + B * sAPTotal[i] );
    }
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
void Shift_Arr()
{
  for (int i = 0; i < nDeme - 1; i++)
    {
      dA[i][0] = dA[i+1][0] ;
      dA[i][1] = dA[i+1][1] ;
    }
  
  dA[nDeme - 1][0] = 0 ;
  dA[nDeme - 1][1] = 0 ;

  pDisp ++ ;
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
void Update_Deme_Det () 
{
  // Deterministic Update of Demes from Specie-Fraction
  for (int i = 0; i < nDeme; i++ )
    {
      dA[i][0] = static_cast<int>(sA[i][0] * populationSize) ;
      dA[i][1] = static_cast<int>(sA[i][1] * populationSize) ;
    }
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
void Update_Deme_DU() 
{
  // Deterministic Update (DU) for total population
  // Stochastic Update for strains
  long double spFrac, p ; 
  long spNum ; 

  for (int i = 0 ; i < nDeme ; i++ ) 
    {
      spFrac = sA[i][0] + sA[i][1] ;
      p = sA[i][0] / spFrac ; 
      spNum = spFrac * populationSize ; 
      
      dA[i][0] = gsl_ran_binomial(r, p, spNum ) ; 
      dA[i][1] = int(spNum) - dA[i][0] ; 
    }
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
void Update_Deme_SU() 
{
  // Fully Stochastic Update (SU) of front
  double pAux[nSp + 1] ;
  unsigned int nAux[nSp + 1] ;

  for (int i = 0 ; i < nDeme ; i++ ) 
    {
      pAux[0] = sA[i][0] ; 
      pAux[1] = sA[i][1] ; 
      pAux[nSp] = 1.0 - pAux[0] - pAux[1] ; 

      gsl_ran_multinomial(r, nSp + 1, populationSize, pAux, nAux ) ;  // Obtain Next Gen. 

      dA[i][0] = nAux[0] ; 
      dA[i][1] = nAux[1] ; 
    }
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
double calcFGlobal ( long arr[][nSp], int n) 
{
  // Calculates the fraction of first strain across the front
  // Input Specie Array and nDeme
  double f = 0.0 ; 
  double population = 0.0 ; 
  for (int i = 0; i < n; i++ ) 
    {
      population += double( arr[i][0] + arr[i][1] ) ; 
      f += double( arr[i][0] ) ; 
    }

  return f/population ; 
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
double calcHetAverage ( long arr[][nSp], int n) 
{
  // Calculates global heterozygosity using total number of each strain across the front
  // Input Specie Array and nDeme
  double h = 0.0 ; 
  int counter = 0 ; 
  double population ; 

  for (int i = 0 ; i < n ; i++ ) 
    {
      population = double( arr[i][0] ) + double( arr[i][1] )  ;
      
      if (population > 0.0) 
	{
	  h += ( arr[i][0] / population ) * (1.0 - arr[i][0] / population ) ; 
	  h += ( arr[i][1] / population ) * (1.0 - arr[i][1] / population ) ; 
	  
	  counter ++ ; 
	}	
    }
  return h / counter ; 
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
double calcAverageHet (long arr[][nSp], int n ) 
{
  // Calculates average heterozigosity by averaging local values of H across the front
  // different quantity than calcHetAverage
  double population = 0.0 ; 
  double p = 0.0 ; 
  for (int i = 0; i < n; i++ ) 
    {
      population += float( arr[i][0] + arr[i][1] ) ; 
      p += arr[i][0] ; 
    }
  return 2.0 * p/population * (1 - p/population) ; 
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
void HetZ_Update() 
// Adds current values to HetZ vector 
{
  hetGlobal.push_back( calcHetAverage( dA, nDeme )) ; 
  fGlobal.push_back( calcFGlobal( dA, nDeme )) ; 
  pop.push_back( instDeme + pDisp  ) ; 
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
int checkFixation(long arr[][nSp]) 
// Records the result of the fixation experiment 
{
    int result = 0 ; 
    int occupiedDemes = 0 ;
    int speciesNumber = 0 ;
    for (int i = 0; i < nDeme; i++)
    {
        if (arr[i][0] + arr[i][1] != 0)
            occupiedDemes ++ ;
        if (arr[i][0] == 0 and arr[i][1] != 0)
            speciesNumber += 1 ;
        else if (arr[i][0] != 0 and arr[i][1] == 0)
            speciesNumber -= 1 ;
    }

    if (speciesNumber == -occupiedDemes){
        result = 0 ;
    }
    else if (speciesNumber == occupiedDemes){
        result = 1 ;
    }
    else{
        // No fixation
        result = -1 ;
    }

    return result ;
}
// ------------------------------------------------------------------- //
//
//
void simulationLoop(int t)
  {
  Mig_del_1();
  Growth(); 
  if (deterministicUpdate == 0 )
    Update_Deme_SU();
  if (deterministicUpdate == 1 ) 
    Update_Deme_DU();

  instDeme = sumArr( dA, populationSize );     
  while ( instDeme  > nDeme / 2.0 ) 
    {
    Shift_Arr();
    instDeme = sumArr( dA, populationSize ); 
    }
    
  if (t % tHetInterval == 0 ) 
    HetZ_Update();
  }


void runFixedTime(string termination)
  {
  calcNGen(2.0, 4.0) ; // Calculate nGen 
  Initialize();
      
  for (unsigned int t = 0; t < nGen ; t++ )
    {   
    simulationLoop(t) ;
      // ~~~ Profile ~~~ // 
    if (t % tProfInterval == 0) 
      {
      ostringstream profName;
      profName << "profile_snapshots_N" << populationSize << "_r" << r0 << "_m" << m0 << "_A" << A << "_B" << B << 
          "_run" << runNumber << termination ; 
      ofstream fprofTime;
      fprofTime.open( profName.str(), std::fstream::app ) ; 

      fprofTime << t ;
      for ( int k = 0; k < nDeme; k++ )
        fprofTime << ',' << dA[k][0] << ',' << dA[k][1] ;
      fprofTime << endl ; 
      fprofTime.close() ;
      }
     // ~~~~~~~~~~~~~~ // 
    }
  }


// ################################################################## // 
// ################################################################## // 


// ------------------------------------------------------------------- //
int main(int argc, char * argv[] ) 
  {

  // ~~~ RNG Env Setup ~~~ // 
  gsl_rng_env_setup();
  T = gsl_rng_mt19937 ;
  r = gsl_rng_alloc(T) ; 
  // ~~~~~~~~~~~~~~~~~~~ // 

  // ~~~ SEED Set-Up ~~~ //
  int sysRandom ; 
  ifstream firand ; 
  firand.open("/dev/urandom", ios_base :: in | ios_base :: binary) ; 
  firand.read((char*) (&sysRandom), sizeof sysRandom ) ; 
  firand.close();
    cout << "Random seed:" << sysRandom << endl; 
  gsl_rng_set(r, sysRandom) ; 
  // ~~~~~~~~~~~~~~~~~~ // 


  // ~~~ Input Parameters ~~~ // // UPDATE for A. 
  int c ; 
  bool fix = false ; // Run fixation simulation
  runNumber = 0 ; 

  while ( ( c = getopt ( argc, argv, "n:r:B:N:m:A:d:w:x:b:" ) ) != -1 ) // Params //
    {
      // 'm'> 'm0'  ; p >> A. Since only single quotes can be used in follow sections
      // -xx IMPT
      if (c == 'n' ) 
	    runNumber = atoi (optarg) ; 
      else if ( c == 'r' ) 
	    r0 = atof(optarg) ; 
      else if ( c == 'B' ) 
	    B = atof(optarg) ; 
      else if ( c == 'N' )
	    populationSize = atof(optarg) ; 
      else if ( c == 'm' )
	    m0 = atof(optarg) ;
      else if ( c == 'A' ) 
	    A = atof(optarg) ; 
      else if ( c == 'd' ) 
	    deterministicUpdate = atof(optarg) ;
      else if ( c == 'w' ) 
        // Not implemented yet 
	    widthCheckFlag = atof(optarg) ;
      else if ( c == 'x' )
        {
        xFixation = atoi(optarg) ;
        fix = true ;
        }
      else if (c == 'b' )
        batchNumber = atoi( optarg ) ;
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~ //  
  // Check values are consistent
  double maxMigrationRate = m0*( 1.0 + A ) ;
  // Set maximum migration rate to 0.3
  if ( maxMigrationRate > 0.3 )
    {
    A = 0.3/m0 - 1.0;
    cout << "Error! Maximum migration rate is too high. A set to " << A << endl ;
    }
      
  // ~~~ Out Files Name Variables ~~~ //
  ostringstream str_populationSize ; 
  ostringstream str_r0 ; 
  ostringstream str_m0 ; 
  ostringstream str_A ; 
  ostringstream str_B ; 
  ostringstream str_run ; 
  ostringstream str_nDeme ; 

  // Keep trailing zeros for integer values
  double intPart ;
  double fractPart ;
  fractPart = modf( A , &intPart ) ;
  if( !fractPart )
    {
    str_A << fixed << setprecision(1) ; 
    }

  fractPart = modf( B , &intPart ) ;
  if( !fractPart )
    {
    str_B << fixed << setprecision(1) ;
    }

  // Assign values to streams
  str_populationSize << populationSize ;
  str_r0 <<  r0; 
  str_m0 << m0 ; 
  str_A << A ; 
  str_B << B ; 
  str_run << runNumber ; 
  str_nDeme << nDeme ; 
  string termination = "_demes" + str_nDeme.str() + ".txt" ;
  // ~~~~~~~~~~~~~~~~~~~~~//


  // ~~~ Simulation ~~~ // 
  clock_t c_init = clock() ; 
  if( fix == true )
    {
    cout << "Starting fixation simulations..." << endl ;

    calcNGen(2.0, 4.0) ;

    // Check for input file and initialize
    ostringstream finput_name ;
    finput_name << "profile_N" << populationSize << "_r" << r0 << "_m" << m0 << "_A" << str_A.str() << "_B" << str_B.str() << "_avg.txt" ;
    int initializeError = fixationInitialize( xFixation, finput_name.str() ) ;
    if ( initializeError )
        {
            cout << "ERROR! Front initialization failed. Quitting program" << endl ;
            return -2 ;
        }

    int t = 0 ;
    while( checkFixation( dA ) == -1 )
      {   
      simulationLoop(t) ;
      t++ ;
      }

    int fixedStrain = checkFixation( dA ) ;
    cout << "Fixed strain: " << fixedStrain << endl;

    // Output fixation result
    ostringstream str_x ; 
    ostringstream str_batch ; 

    str_x << xFixation ;
    str_batch << batchNumber ;

    string fixName = "fixation_N" + str_populationSize.str() + "_r" + str_r0.str() + "_m" + str_m0.str() + "_A" + str_A.str() + "_B" + str_B.str() + "_x" + str_x.str() + "_batch" + str_batch.str() + termination ;
    ofstream ffix ;
    ffix.open( fixName, fstream::app ) ;
    ffix << runNumber << "," << xFixation << "," << fixedStrain << endl ;
    ffix.close() ;
    }
  else
    {
    cout << "Starting fixed time simulations w/ nGen = " << nGen << endl ; // ~~~ nGen is wrong at this point
    runFixedTime( termination ) ;
    }
  clock_t c_fin = clock() ; 
  // ~~~~~~~~~~~~~~~~~~~~~~~ //


  // ~~~ OUTPUT ~~~ //
  string velName = "velocity_N" + str_populationSize.str() + "_r" + str_r0.str() + "_m" + str_m0.str() + "_A" + str_A.str() + "_B" + str_B.str() + "_run" + str_run.str() + termination ; 
  string hetName = "hetero_N" + str_populationSize.str() + "_r" + str_r0.str() + "_m" + str_m0.str() + "_A" + str_A.str() + "_B" + str_B.str() + "_run" + str_run.str() + termination ; 
  string profName = "profile_N" + str_populationSize.str() + "_r" + str_r0.str() + "_m" + str_m0.str() + "_A" + str_A.str() + "_B" + str_B.str() + "_run" + str_run.str() + termination ; 
  string paramName = "param_N" + str_populationSize.str() + "_r" + str_r0.str() + "_m" + str_m0.str() + "_A" + str_A.str() + "_B" + str_B.str() + termination ; 

  ofstream fpop ;
  ofstream fhet ;
  ofstream fprof ;
  ofstream fparam ; 

  fhet.open(hetName) ; 
  fpop.open(velName) ; 
  fprof.open(profName) ; 
  fparam.open(paramName) ; 
  // ~~~~~~~~~~~~~~~~~~~~~~~~ //


  fparam << "Number of generation:" << nGen << ", Number of species:"<< nSp << ", Population Size:" << populationSize << ", Growth rate:" << r0 << ", Growth cooperativity:" << B << ", Migration rate:" << m0 << ", Migration cooperativity:" << A << ", Number of demes:" << nDeme << ", NoiseFlag:" << deterministicUpdate  ; 
  fparam.close() ;

  // ~~~ Hetrozygocity
  for (unsigned int i = 0; i < hetGlobal.size(); i++ ) 
    fhet << tHetInterval * i << "," << hetGlobal[i] << "," << fGlobal[i] << endl; 
  fhet.close() ;

  // ~~~ Population
  for (unsigned int i = 0; i < pop.size(); i++ ) 
    fpop << tHetInterval * i << "," << pop[i] << endl; 
  fpop.close() ; 

  // Profile 
  for ( int i = 0 ; i < nDeme ; i++ ) 
	  fprof << i << "," << dA[i][0] << "," << dA[i][1] << endl ; 
  fprof.close() ;

  // ~~~~~~~~~~~~~ //

  cout << "Simulation time: " << nGen << endl;
  cout << "Initial Population: " << pop.front() << '\t' << "Final Population: " << pop.back() << endl ; 
  cout << "Estimated velocity: " << ( pop.back() - pop.front() )/nGen << endl ; // Not accurate but sufficient 
  cout << endl << "Final heterozygosity: " << calcHetAverage( dA, nDeme) << endl ; 
  cout << "Run time: " << double(c_fin - c_init ) / CLOCKS_PER_SEC << endl ;  

  return 0;
  } // END MAIN

// ------------------------------------------------------------------- //
