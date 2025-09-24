#include "Genetic_algorithm.h"
#include "random.h"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <mpi.h>
using namespace std;
using namespace arma;

int main(int argc, char *argv[]) {
//initialize the parallel process
  MPI_Init( &argc, &argv);

//store rank and number of processes
int rank, size;

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

//generate a population for each continent (process)

  int N_cities = 110;
  int N_chromosomes=200;
  int N_migr = 100;
  int Ngen=2000;
  int Ntot=Ngen/N_migr;
  int migrants=4; //15
  Random rnd;
  int bestpath[N_cities];

   int seed[4];
  int p1, p2;
  ifstream Primes("../INPUT/Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("../INPUT/seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  rnd.SaveSeed();

  string filename="../INPUT/cap_prov_ita.dat";
  GA Darwin(N_chromosomes, N_cities, rank);
  mat cities= Load_cities(filename, N_cities);
  //mat cities=Generate_cities(N_cities, 0);
  mat D=Distances(cities);
  ofstream out;
  ofstream outbest;


  outbest.open("../OUTPUT/Best_loss"+to_string(rank)+".dat");
    
  out.open("../OUTPUT/besthalfloss"+to_string(rank)+".dat");
   out<<"#gen      BHL"<<endl;
  

for(int k=0;k<Ntot;k++){
  if(rank==0)
  cout<<"migration "<<k+1<<endl;

 //perform indepent GA searches until migration
  for(int i=0; i<N_migr;i++){
    
    Darwin.Evolve(D);
    
      Darwin.Sort_by_fitness(D);
     out<<setw(12)<<i<<setw(12)<<Darwin.BestHalfLoss(D)<<endl;
    outbest<<setw(12)<<i<<setw(12)<<Darwin.Loss(D,0)<<endl;
    
  }
  //every continent performs a migration of migrants individuals (which are the best migrants one of each population)
  for(int l=0;l<migrants;l++){
 //store the best path in an int array
  Darwin.Sort_by_fitness(D);
  for(int i=0;i<N_cities;i++){
     bestpath[i]=Darwin.Get_Chromosome(l).Get_Gene(i);
  }
  
  //the first process chooses a random permutation of the ranks and communicate it to all of the continents
    int migration[size];
  if(rank==0){

  for(int i=0;i<size;i++){
    //fill the array with the positions
    migration[i]=i;
  }
  //shuffle the array
  int temp=0;
    for(int i=size-1;i>0;i--){
      int j= (int)rnd.Rannyu(0.0,(double)i);
    temp=migration[i];
    migration[i]=migration[j];
    migration[j]=temp;
  }
  //communicate the protocol to everyone
  }
   MPI_Bcast(migration, size, MPI_INT, 0, MPI_COMM_WORLD);

  //now everyone in the list performs the migration by exchanging its population with his neighbor in migration array
  int partner=-1; //to avoid selecting automatically the first one
for(int i=0; i<size-1;i+=2){
  //everyone stores its partner
  if(migration[i]==rank)
  partner=migration[i+1];
  //symmetric case
  if(migration[i+1]==rank)
  partner=migration[i];
}
//prepare the array to receive the population
  int recvpath[N_cities];
  if(partner != -1){ //if the number of processes isn't even, one of them simply doesn't migrate

// set a condition to make sure that one sends first
  if(rank<partner){ //se passo un array non serve & perche è gia un puntatore
    MPI_Send(bestpath, N_cities, MPI_INT, partner, 0, MPI_COMM_WORLD);
    MPI_Recv(recvpath, N_cities, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Darwin.Replace(recvpath, D);
    //cout << "Loss after migration on rank " << rank <<" to "<<partner<< ": " << Darwin.Loss(D, 0) << endl;

    
  }
  else{
    MPI_Recv(recvpath, N_cities, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(bestpath, N_cities, MPI_INT, partner, 0, MPI_COMM_WORLD);
    Darwin.Replace(recvpath, D);
  }

  } 
  }

}
//after migrations, the first process gathers all of the best paths and chooses the best one
 //store the best path in an int array
  Darwin.Sort_by_fitness(D);
  for(int i=0;i<N_cities;i++){
     bestpath[i]=Darwin.Get_Chromosome(0).Get_Gene(i);
  }
 int bestpaths[size*N_cities];
//only the first process must have the gathering variable
  MPI_Gather(bestpath, N_cities, MPI_INT, bestpaths, N_cities, MPI_INT, 0, MPI_COMM_WORLD);
if(rank==0){
  for(int i=0;i<size;i++){
    int best[N_cities];
    for(int j=0;j<N_cities;j++){
      best[j]=bestpaths[i*N_cities+j];
     
  }
  Darwin.Replace(best, D);

}
  Darwin.BestHalfLoss(D);
  Darwin.Pick_best(D);

}

  outbest.close();
  out.close();

MPI_Finalize();
  return 0;

}