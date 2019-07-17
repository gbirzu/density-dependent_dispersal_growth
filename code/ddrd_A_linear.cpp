// Code Version 5.x C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
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

using namespace std;

// ~~~ Variables ~~~ //
const double pi = 3.14159 ;
const int nDeme = 300 ;
const int nSp = 2 ;
unsigned int nGen = 1000 ; 

double gR = 0.01 ;
double m0 = 0.1 ;
double m1 = 0.0 ;
int Kap = 1000 ;  // Carrying Capacity
int pDisp, N ;
int run ;          // Run Number 
double delta ;     // Raises to non-linear power
double fStar = 0.0 ; // Alee Effect (UNUSED) 

int nHetPoints = 1000 ; 
int tHetInterval = int( nGen / nHetPoints) ; 
int widthCounter = 0 ;  // 
// ~~~~~~~~~~~~~~~~ //

// ~~~ Flag Variables ~~~ //
int instDeme , maxDeme ; 
unsigned int noiseFlag = 1 ; // 0: _DF || 1: _EN
unsigned int checkLastDeme = 1; 
unsigned int widthCheckFlag = 0 ; 
// ~~~~~~~~~~~~~~~~~~~~~ //


// ~~~ Array ~~~ //
long double sA[nDeme][nSp], sAP[nDeme][nSp], sAPTotal[nDeme] ;
long int dA[nDeme][nSp] ;

vector<double> pop ; 
vector<double> het ; 
vector<double> hetOld ; 
// ~~~~~~~~~~~~~ //


// ~~~ RND GSL ~~~ //
const gsl_rng_type * T ; 
gsl_rng * r ; 
// ~~~~~~~~~~~~~~ //



// ~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~ //

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
      dA[i][0] = gsl_ran_binomial(r, 0.5, Kap) ;
      dA[i][1] = Kap - dA[i][0] ;

      sA[i][0] = dA[i][0] / float(Kap) ;
      sA[i][1] = dA[i][1] / float(Kap) ;
    }

  pDisp = 0 ; 
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
void calc_nGen ()
{
  double t_bsc ;
  t_bsc = (pow(log(Kap),3))/(2*gR*pow(pi, 2)) ;
  double mRatio ;
  if (m0 == 0)
      mRatio = 10 ;
  else
      mRatio = m1/m0 ;
  if (mRatio < 1.)
      nGen = max(10*Kap, max(int(10*t_bsc), 10000));
  else if (mRatio >= 1. and mRatio < 2.)
  {
      double exponent = mRatio - 1. ;
      nGen = max(10000, max(int(t_bsc), 10*int(sqrt(0.25/(2*gR))*pow(Kap, exponent)))) ;
  }
  else
      nGen = 20*Kap ;

  tHetInterval = int( nGen / nHetPoints) ; 
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
double sumArr( long  arr[][nSp], int Kap) // Only for dA
{
  double sum = 0.0 ; 
  for (int i = 0 ; i < nDeme ; i++)
    {
      sum += arr[i][0] / double(Kap) + arr[i][1] / double(Kap) ; 
    }
  return sum; 
}
// ------------------------------------------------------------------ //



// ------------------------------------------------------------------- //
double mig_rate( double sp_frac) // defines function for migration rate; uses m0 and m1 as parameters
{
    return m0 + m1*sp_frac;
}
// ------------------------------------------------------------------- //


// ------------------------------------------------------------------- //
void Mig_del_1() // D_1 term is polynomial of order 1. For Faster simulations 
{

  for (int i = 0; i < nDeme; i++) 
    {
      // Converts into Fractions
      sA[i][0] = double(dA[i][0]) / double(Kap) ;
      sA[i][1] = double(dA[i][1]) / double(Kap) ;

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
void Growth_Logistic()
{
  for (int i = 0; i<nDeme; i++) 
    {
      // Assign Prev. Array
      sAP[i][0] = sA[i][0] ;
      sAP[i][1] = sA[i][1] ;
      sAPTotal[i] = sAP[i][0] + sAP[i][1] ;
      

      // Implement Growth; per capita growth depends only on total number since species are neutral
      sA[i][0] = sAP[i][0] + gR * ( 1.0 - sAPTotal[i] ) * sAP[i][0] ;
      sA[i][1] = sAP[i][1] + gR * ( 1.0 - sAPTotal[i] ) * sAP[i][1] ;
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
void Update_Deme_Det () // Deterministic Update of Demes from Specie-Fraction
{
  for (int i = 0; i < nDeme; i++ )
    {
      dA[i][0] = static_cast<int>(sA[i][0] * Kap) ;
      dA[i][1] = static_cast<int>(sA[i][1] * Kap) ;
    }
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
void Update_Deme_DF() 
// Only Demographic Fluctuations. NO Ecological Noise
{
  long double spFrac, p ; 
  long spNum ; 

  for (int i = 0 ; i < nDeme ; i++ ) 
    {
      spFrac = sA[i][0] + sA[i][1] ;
      p = sA[i][0] / spFrac ; 
      spNum = spFrac * Kap ; 
      
      dA[i][0] = gsl_ran_binomial(r, p, spNum ) ; 
      dA[i][1] = int(spNum) - dA[i][0] ; 
    }
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
void Update_Deme_EN() 
// Multinomial sampling (Ecological Noise) for genetic drift
{
  double pAux[nSp + 1] ;
  unsigned int nAux[nSp + 1] ;

  for (int i = 0 ; i < nDeme ; i++ ) 
    {
      pAux[0] = sA[i][0] ; 
      pAux[1] = sA[i][1] ; 
      pAux[nSp] = 1.0 - pAux[0] - pAux[1] ; 

      gsl_ran_multinomial(r, nSp + 1, Kap, pAux, nAux ) ;  // Obtain Next Gen. 

      dA[i][0] = nAux[0] ; 
      dA[i][1] = nAux[1] ; 
    }
}
// ------------------------------------------------------------------- //



// ------------------------------------------------------------------- //
double calcHet ( long arr[][nSp], int n) 
// Input Specie Array and nDeme. 
{
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
double calcHetOld (long arr[][nSp], int n ) 
{
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
  het.push_back( calcHet( dA, nDeme )) ; 
  hetOld.push_back( calcHetOld( dA, nDeme )) ; 
  pop.push_back( instDeme + pDisp  ) ; 
}
// ------------------------------------------------------------------- //




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


  // ~~~ Input Parameters ~~~ // // UPDATE for m1. 
  int c, run = 0 ; 
  while ( (c = getopt (argc, argv, "r:f:m:p:g:K:n:w")) != -1 ) // Params //
    // 'm'> 'm0'  ; p >> m1. Since only single quotes can be used in follow sections
    // -xx IMPT

    {
      if (c == 'r' ) 
	run = atoi (optarg) ; 
      //else if (c == 'f' ) 
      //fStar = atof( optarg) ; 
      else if ( c == 'm' )
	m0 = atof(optarg) ;
      else if ( c == 'p' ) 
	m1 = atof(optarg) ; 
      else if ( c == 'g' ) 
	gR = atof(optarg) ; 
      else if ( c == 'K' )
	Kap = atof(optarg) ; 
      else if ( c == 'n' ) 
	noiseFlag = atof(optarg) ;
      else if ( c == 'w' ) 
	noiseFlag = atof(optarg) ;
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~ //  


  // ~~~ Out File Name Var ~~~ //
  ostringstream str_Kap ; 
  str_Kap << Kap ;

  ostringstream str_gR ; 
  str_gR <<  gR; 
  
  ostringstream str_m0 ; 
  str_m0 << m0 ; 

  ostringstream str_m1 ; 
  str_m1 << m1 ; 

  ostringstream str_run ; 
  str_run << run ; 

  ostringstream str_nDeme ; 
  str_nDeme << nDeme ; 

  string termination = "_demes" + str_nDeme.str() + ".txt" ;

  string velName = "velocity_N" + str_Kap.str() + "_gf" + str_gR.str() + "_m0" + str_m0.str() + "_m1" + str_m1.str() + "_run" + str_run.str() + termination ; 
  string hetName = "hetero_N" + str_Kap.str() + "_gf" + str_gR.str() + "_m0" + str_m0.str() + "_m1" + str_m1.str() + "_run" + str_run.str() + termination ; 
  string profName = "profile_N" + str_Kap.str() + "_gf" + str_gR.str() + "_m0" + str_m0.str() + "_m1" + str_m1.str() + "_run" + str_run.str() + termination ; 
  string paramName = "param_N" + str_Kap.str() + "_gf" + str_gR.str() + "_m0" + str_m0.str() + "_m1" + str_m1.str() + "_run" + str_run.str() + termination ; 

  ofstream fpop, fhet, fprof, fparam ; 

  fparam.open(paramName) ; 
  fhet.open(hetName) ; 
  fpop.open(velName) ; 
  fprof.open(profName) ; 
  // ~~~~~~~~~~~~~~~~~~~~~~~~ //


  // ~~~ Simulation ~~~ // 
  clock_t c_init = clock() ; 
  int count ; 
  
  calc_nGen() ; // Calculate nGen 
  Initialize();
      
  for (unsigned int t = 0; t < nGen ; t++ )
    {   
      Mig_del_1() ;
      Growth_Logistic() ; 
      if (noiseFlag == 0 )
	{ Update_Deme_DF() ; }
      if (noiseFlag == 1 ) 
	{ Update_Deme_EN() ; }

     instDeme = sumArr( dA, Kap ) ;     
     while ( instDeme  > nDeme / 2.0 ) 
       {
	 Shift_Arr();
	 instDeme = sumArr( dA, Kap ) ; 
       }
     count = t ; 
    
     if (t % tHetInterval == 0 ) 
       {HetZ_Update();  }

     // ~~~ Profile ~~~ // 
     if ( widthCheckFlag != 0 ) 
       if (t == pow(2, widthCounter)) 
	 {
	   ostringstream str_t ; 
	   str_t << t ; 
	   string profName = "profile_t" + str_t.str() + "_N" + str_Kap.str() + "_gf" + str_gR.str() + "_mig0" + str_m0.str() + "_mig1" + str_m1.str() + "_run" + str_run.str() + termination ; 
	   ofstream fprofw ; // Terrible file name
	   fprofw.open( profName ) ; 

	   for ( int k = 0; k < nDeme; k++ )
           fprofw << k << ',' << dA[k][0] << ',' << dA[k][1] << endl ; 
	   widthCounter++ ; 

	 }
     // ~~~~~~~~~~~~~~ // 

    }
  clock_t c_fin = clock() ; 
  // ~~~~~~~~~~~~~~~~~~~~~~~ //


  // ~~~ OUTPUT ~~~ //

  fparam << "Number of generation:" << nGen << ", Number of species:"<< nSp << ", Population Size:" << Kap << ", Growth rate:" << gR << ", Migration_0 rate:" << m0 << ", Migration_1 rate:" << m1 << ", Number of demes:" << nDeme << ", NoiseFlag:" << noiseFlag  ; 

  // ~~~ Hetrozygocity
  for (unsigned int i = 0; i < het.size(); i++ ) 
    fhet << tHetInterval * i << "," << het[i] << "," << hetOld[i] << endl; 

  // ~~~ Population
  for (unsigned int i = 0; i < pop.size(); i++ ) 
    fpop << tHetInterval * i << "," << pop[i] << endl; 
  fpop.close(); 

  // Profile 
  for ( int i = 0 ; i < nDeme ; i++ ) 
	  fprof << i << "," << dA[i][0] << "," << dA[i][1] << endl ; 

  // ~~~~~~~~~~~~~ //

  cout << "Simulation time: " << nGen << endl;
  cout << "Initial Population: " << pop.front() << '\t' << "Final Population: " << pop.back() << endl ; 
  cout << endl << "Final heterozygosity: " << calcHet( dA, nDeme) << endl ; 
  cout << "Run time" << double(c_fin - c_init ) / CLOCKS_PER_SEC << endl ;  

  return 0;
} // END MAIN

// ------------------------------------------------------------------- //
