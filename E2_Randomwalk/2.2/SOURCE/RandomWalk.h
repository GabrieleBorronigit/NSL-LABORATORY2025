#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include "random.h"
#include <string>
#include <vector>
using namespace std;


class RandomWalk{

private:

int m_D;                                //number of dimensions
int m_N_steps;                          //number of steps
double m_a;                             //step length
vector<double> m_r;                     //position vector
vector<vector<double>> m_History;       //matrix that saves every step
Random m_gen;                          //RNG 

vector<double>  m_amsd_lattice;
vector<double> m_amsd_cont;



public:


RandomWalk(int N_steps, double a, Random rnd, int D=3): //constructor, by default sets number of dimensions to three and the starting position in the origin
 m_D{D}, m_N_steps{N_steps}, m_a{a}, m_r(D,0), m_gen{rnd}{}



  void Set_Steps(int N_steps){ m_N_steps=N_steps;}          //sets the number of steps for the walk

  void Set_Step_Length(double a){m_a=a;}                    //sets the step length


  void Reset(){m_r=vector<double>(m_D,0);                   //resets the walk by putting the walker in the origin
  m_History.clear();
  m_amsd_lattice.clear();
  m_amsd_cont.clear();

  }


   void PrintPosition() const {                             //prints current position (debug)
        for (int i = 0; i < m_D; ++i) cout << m_r[i] << " "<<endl;
    }


double Get_amsd_lattice_comp(int step) const{
    return m_amsd_lattice[step];
}


double Get_amsd_cont_comp(int i){
    return m_amsd_cont[i];
}

void History_ToFile(string filename);                       //prints the walk history



double Distance();                                          //evaluates the distance from the origin






//methods for walking

void LatticeWalk(bool save_history=false); //saving history is optional, since it wasn't explicitly requested, it doesn't do so by default





vector<double> Uniform_continuum_direction(); //Samples random direction in the continuum

void ContinuumWalk(bool save_history=false); //self explainatory



void Save_History(){m_History.push_back(m_r);} //saves current position in history



};



//data blocking methods

double error(double av, double av2, int n);



void ComputeMeanErrors_DB(const vector<double> &data, int N_blocks, ofstream &out, int step);