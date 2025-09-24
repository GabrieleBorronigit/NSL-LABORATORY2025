#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random.h"
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm> //for min()
#include "VariationalMC.h"

using namespace std;



void VMC::initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  ofstream couth("../OUTPUT/local_energy.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couth << "#   N_BLOCK:  <H>Block    block ave     error:" << endl;
  couth.close();

  ofstream coutpdf("../OUTPUT/PDF.dat");
  coutpdf << "#N_BLOCK:  bin     ave     prog ave     error:" << endl;
  couth.close();
  //initialize arbitrary mu and sigma and random starting position
  _mu = 1.0;
  _sigma = 0.6;
  _x = _mu + _sigma * _rnd.Rannyu(-1,1);  // Start inside the wavefunction support

  //initialize acceptance
  _naccepted=0;
  _nattempts=0;
  //initialize data blocking variables
  _block_av = 0;
  _average = 0;
  _global_av = 0;
  _global_av2 = 0;

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("OUTPUT/output.dat");
  string property;
  
  while ( !input.eof() ){
    input >> property;
   
       coutf << "VARIATIONAL MONTE CARLO (VMC) SIMULATION"  << endl;
 
     if( property == "DELTA" ){
      input >> _delta;
      coutf << "DELTA= " << _delta << endl;}
    else if( property == "DELTAMU" ){
      input >> _delta_mu;
      coutf << "DELTAMU= " << _delta_mu << endl;}
    else if( property == "DELTASIGMA" ){
      input >> _delta_sigma;
      coutf << "DELTASIGMA= " << _delta_sigma << endl;}
    else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "TEMPERATURE"){
        input>>_T;
        _beta=1.0/_T;
        coutf<< "STARTING TEMPERATURE="<< _T <<endl;
    } else if( property == "BINS"){
      input>>_nbins;
      xmin=-3.0;
      double xmax=3.0;
      _bin_size=(xmax-xmin)/_nbins;
      pdf.resize(_nbins);

    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else {cerr << "PROBLEM: unknown input" << endl;
    exit(1);
    return;}
  }
  input.close();

  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}






double VMC:: Psi_T() {
    return exp(-pow((_x - _mu), 2) / (2 * _sigma * _sigma)) +
           exp(-pow((_x + _mu), 2) / (2 * _sigma * _sigma));
}

double VMC:: Psi_T(double x) {
    return exp(-pow((x - _mu), 2) / (2 * _sigma * _sigma)) +
           exp(-pow((x + _mu), 2) / (2 * _sigma * _sigma));
}

double VMC::Psi_2(){
return Psi_T()*Psi_T();

}
double VMC::Psi_2(double x){
return Psi_T(x)*Psi_T(x);

}


double VMC::V(){
    return pow(_x,4)-5.0/2.0*pow(_x,2);
}

void VMC::Move(){
 //proposed move
    double xnew= _x+ _delta*_rnd.Rannyu(-1.0,1.0);
    //evaluate Metropolis acceptance
    _nattempts++;
    if(Metropolis(_x, xnew))
    _x=xnew;
}

bool VMC::Metropolis( double xold, double xnew){
    bool decision= false;
    //Metropolis acceptation step
    double acceptance= min(1.0, Psi_2(xnew)/Psi_2(xold));
    if(_rnd.Rannyu()<=acceptance){
    decision=true;
    _naccepted++;}
   
    return decision;
}


double VMC::Kinetic_Energy(){
return 
-0.5*(exp(-pow((_x - _mu), 2) / (2.0 * pow(_sigma, 2))) * 
    (pow((_x - _mu), 2) / pow(_sigma, 4) - 1.0 / pow(_sigma, 2)) +
exp(-pow((_x + _mu), 2) / (2.0 * pow(_sigma, 2))) * 
    (pow((_x + _mu), 2) / pow(_sigma, 4) - 1.0 / pow(_sigma, 2)))/Psi_T();
}

void VMC::Trial_Mean_Energy(){
//one sample of <H>
double H=0;
for(int i=0; i<_nsteps; i++){
//sample ddp psi^2
Move();
H+=(Kinetic_Energy()+V());
}
_block_av+=H/_nsteps;
return;
}


pair<double, double> VMC::Trial_Mean_Energy_SA(){ //returns [<H>, error]
//one sample of <H>
double H=0, H2=0;
for(int i=0; i<_nsteps; i++){
//sample ddp psi^2
Move();
H+=(Kinetic_Energy()+V());
H2+=pow((Kinetic_Energy()+V()),2);
}
return make_pair(H/_nsteps, this->error(H,H2, _nsteps));
}


void VMC:: PDF(int Nsamples){
  //samples |psi|^2 as an histogram, after optimization
  ofstream out("../OUTPUT/PDF.dat", ios::app);
for(int i=0;i<Nsamples;i++){
  Move();
  int bin=int((_x-xmin)/_bin_size);
  if(bin<_nbins&& bin>=0)
  pdf[bin]+=1;
}

//normalize
for(int i = 0; i < _nbins; i++) {
//normalization in python
  out<<pdf[i]<<endl;
}

out.close();

}
void VMC :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av=0;
  _naccepted=0;
  _nattempts=0;
  _x = _mu + _sigma * _rnd.Rannyu(-1.0,1.0);  // random intorno a ±σ

  return;
}



void VMC :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av;
  _global_av  += _average;
  _global_av2 += _average * _average;

    coutf.open("../OUTPUT/local_energy.dat",ios::app);
    average  = _average;
    sum_average = _global_av;
    sum_ave2 = _global_av2;
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  



   double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}

double VMC :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}



void VMC::Optimize(int SA_steps){
    ofstream coutf;
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Simulated annealing optimization running... " << endl;
    double mu,sigma;
    double mu_old=_mu;
    double sigma_old=_sigma;
    double sigma_best{},mu_best{};
    double mu_global_best=1e6, sigma_global_best=1e6;
    double beta=_beta;
    double cooling=1.005;
    pair <double,double> E_best=make_pair(1e8,1e8), E_global_best=make_pair(1e8,1e8);
    pair <double, double> E=make_pair(0.0,0.0); //to store best energy and error for each temperature to be printed
    //output: acceptance, parameters, energies
    ofstream out;
    out.open("../OUTPUT/acceptance_sa.dat");
    out<<"# T          A"<<endl;


    ofstream outpath;
    outpath.open("../OUTPUT/parameters_walk.dat");
    outpath<<"# T          sigma           mu"<<endl;


    ofstream outen;
    outen.open("../OUTPUT/energies_SA.dat");
    outen<<"#T           <H>          error         SIGMA         MU "<<endl;
    E=Trial_Mean_Energy_SA();
    outen<<setw(12)<<1.0/beta<<setw(12)<<E.first<<setw(12)<<E.second<<setw(12)<<_sigma<<setw(12)<<_mu<<endl; //prnt starting parameters


    for(int i=0; i<SA_steps; i++){ //cycle over desired SA steps
        int naccepted=0; //acceptance evaluation parameters
        int nattempts=0;
 //loop over the temperature, by making more samples in cooler regions: N(T)=N_0/T
        for(int n=0; n<(int)(SA_steps*beta/100);n++){ 
        nattempts++;
  //evaluate current state's energy (and error)
        pair <double, double> Eold=Trial_Mean_Energy_SA();
  do{
        mu=mu_old+_delta_mu*_rnd.Rannyu(-1.0,1.0); //trial random parameter
        sigma=sigma_old+_delta_sigma*_rnd.Rannyu(-1.0,1.0); //sigma must be positive
  }while(sigma<=0.05 || sigma>=2.5 );
        _mu=mu;
        _sigma=sigma;

        pair <double, double> Enew=Trial_Mean_Energy_SA(); //evaluate trial energy

        //metropolis evaluation step
        if(Enew.first>Eold.first){
            //evaluate acceptance probability
            double p=exp(-beta*(Enew.first-Eold.first));
            //accept next configuration with probability p
            if(_rnd.Rannyu()<=p){
                naccepted++;
                mu_old=mu;
                sigma_old=sigma;    
                E=Enew; //store result (<H>_T)
                
            }
            else{   
                _mu=mu_old;
                _sigma=sigma_old;
                E=Eold; //store result
                } }
        else {
       naccepted++;
         mu_old=mu;
        sigma_old=sigma;
        E=Enew;
        }
        //if Enew<Eold there's nothing to be done: the new parameters have already been accepted

    //check if those are the best parameters for the current temperature
    if(E.first<E_best.first){
      E_best=E;
      sigma_best=_sigma;
      mu_best=_mu;
    }
           
        }
        outpath<<setw(12)<<beta<<setw(12)<<sigma_best<<setw(12)<<fabs(mu_best)<<endl;
           out<<setw(12)<<beta<<setw(12)<<(double)naccepted/(double)nattempts<<endl;
           outen<<setw(12)<<beta<<setw(12)<<E_best.first<<setw(12)<<E_best.second<<setw(12)<<sigma_best<<setw(12)<<fabs(mu_best)<<endl;
           if(i%50==0)
           coutf<<"cooling... step #"<<i<<endl;
    //update temperature and reset best variables
   beta*=cooling;
   if(E_best.first<E_global_best.first){
    E_global_best=E_best;
      sigma_global_best=sigma_best;
      mu_global_best=mu_best;
   }
   
  E_best =make_pair(10.0,10.0);  // come reset
  sigma_best = 0.0;
  mu_best = 0.0; 

    }
    //print best results in the output file
    coutf<<"Optimization completed!"<<endl<<"Best Energy and parameters: E="<<E_global_best.first<<"+-"<<E_global_best.second<<"  sigma="<<sigma_global_best<<"  mu="<<fabs(mu_global_best)<<endl;
    coutf.close();
    out.close();
    outpath.close();
    outen.close();
    //_mu=mu_global_best;
    //_sigma=sigma_global_best;
}

void VMC::Loadpars(string filename){
  ifstream in(filename);
  if(!in){
    cerr<<"impossible to open file "<<filename<<endl;
    exit(1);
  }
  double t, sigma, mu;
  while(in>>t>>sigma>>mu){

  _sigma=sigma;
  _mu=mu;
}
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/







