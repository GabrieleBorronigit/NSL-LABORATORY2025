#pragma once
#include "random.h"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
using namespace std;
using namespace arma;


// initialization of TSP
mat Generate_cities(int N, int type);
mat Distances(const mat &cities);
mat Load_cities(string filename, int Ncities);



class Ant{


public:

Ant(int Ncities, Random &rnd_);
void Walk(const mat & Feromons, const mat &D, double alpha, double beta);
int Next_step();
void Reset();
vector<int> Get_visited();
int Get_visited(int i);

private:
int Ncities_; //total length of the trip
vector<int> visited_;
vector<int> unvisited_;
vector<double> Pij_;
Random rnd_;
};


class Colony{


public:
Colony(int Ncities, int Nants, double rho, double Q);
double Loss(const mat&D, int i);
void Update_feromons(const mat&D);
void Search(const mat&D, double alpha, double beta);
void Sort_by_fitness(const mat&D);
void Swap_in_colony(int i, int j);
void Pick_best(const mat&D);
double Best_loss(const mat&D);
double BestHalfLoss(const mat&D);

    private:
    int Ncities_;                   
    int Nants_;                       // #ants
    vector<Ant> Colony_;           //
    mat Feromons_;                   // matrix τ (feromons)
    double rho_, Q_;  // ACO parameters

    vector<int> best_tour_;
    double best_loss_;

    Random rnd_;


};