#pragma once
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random.h"
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm> //for min()


using namespace std;


class VMC{


public:

void initialize(); //initializes the random number generator, data blocking and all the other parameters 
double Psi_T(); //evaluete trial wavefunction
double Psi_T(double x);
double Psi_2(); //square modulus of the trial wavefunction
double Psi_2(double x);
double V();  //evaluate potential energy
double Kinetic_Energy(); //evaluate the fixed cinetic term
void Trial_Mean_Energy(); //measures <H> and accumulates in _block_av
pair <double, double> Trial_Mean_Energy_SA(); //for metropolis acceptance step in simulated annealing
void Move(); //metropolis RW step
bool Metropolis(double xold, double xnew); //evaluate metroolis acceptation step
void block_reset(int blk); //reset block accumulator
void averages(int blk); //data blocking
double error(double acc, double acc2, int blk); //compute errors
int Get_Nbl(){return _nblocks;}
int Get_Nsteps(){return _nsteps;}
void Optimize(int N); //optimizes variational parameters via simulated annealing algorithm
void PDF(int Nsteps); //samples |psi|^2 in an histogram
void Loadpars(string filename);
private:

double _sigma, _mu, _delta_sigma, _delta_mu;//variational parameters
Random _rnd; //generator
int _naccepted, _nattempts; //acceptance parameters
double _x; // coordinate
double _delta; //sets the length of the Metropolis step
double _block_av, _average, _global_av, _global_av2; //data blocking accumulators
int _nsteps, _nblocks; //block size
double _beta, _T; //SA temperature and boltzmann weight
int _nbins;         //for psi sampling
double _bin_size;
double xmin;
vector<double> pdf;

};
