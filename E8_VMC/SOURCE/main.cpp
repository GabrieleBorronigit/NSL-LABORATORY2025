#include "VariationalMC.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random.h"
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm> //for min()
using namespace std;

int main(int argc, char *argv[]){

int type=stoi(argv[1]); //0: optimize + run; 1: only run with best parameters from last run
int opt_steps=1300;

VMC sys;
sys.initialize();           //initialize system's properties (rng, variational parameters, input parameters)
if(type==0){
sys.Optimize(opt_steps);    //run simulated annealing algorithm
sys.block_reset(0);         //reset block accumulators
}
else 
   sys.Loadpars("../OUTPUT/parameters_walk.dat"); 
for(int i=0; i<sys.Get_Nbl();i++){ //measure energies after optimization
    sys.Trial_Mean_Energy();
    sys.averages(i+1);
   sys.block_reset(i+1);
 }
 //sample pdf
 int histSamples=1000000;
 sys.PDF(histSamples);

    return 0;
}