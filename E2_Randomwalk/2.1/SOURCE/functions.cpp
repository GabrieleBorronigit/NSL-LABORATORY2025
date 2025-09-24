#include "functions.h"

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;





double Integral_MC::Uniform_MC(const Template_function &f) {
  //integral estimator
  double G_N{};
  // evaluate the integral
  for (int i = 0; i < m_N; i++) {
    // sample a uniform number in the integration interval
    double x_i = m_a + (m_b - m_a) * m_gen.Rannyu();
    G_N += f.Eval(x_i);
  }
  return (m_b - m_a) * G_N / m_N;
}



double Integral_MC::Importance_Sampling(const Template_function &f, PDF &p) {
  // estimator
  double G_N{};
  // evaluate the integral
  for (int i = 0; i < m_N; i++) {
    // Sample random numbers following given p(x) in the integration interval
    double x_i = p.Sample(m_gen);
    G_N += f.Eval(x_i) / p.Eval(x_i);
  }
  return G_N / m_N;
}






double error(double av, double av2, int n) {  // mean dev std

  if (n == 0)
    return 0;  // set first block error to zero
   

  return sqrt((av2 - pow(av, 2)) / n);
}



//computes statistical averages and errors using data blocking, outputs on file
void ComputeMeanErrors_DB(const vector<double> &data, int N_blocks, string filename) {

  int N_throws = data.size();

  int Blocks_size = static_cast<int>(N_throws / N_blocks);

  double temp_sum{};
  vector<double> average;
  vector<double> square_average;

  //loop over blocks
  for (int i{}; i < N_blocks; i++) {
    temp_sum = 0;
  //loop inside each block
      for (int j{}; j < Blocks_size; j++) {
        temp_sum += data[i * Blocks_size + j];
       }

    average.push_back(temp_sum / Blocks_size);    // average
    square_average.push_back(pow(average[i], 2)); //square average
  }

 //accumulators
  vector<double> mean_prog(N_blocks,0.0);
  vector<double> sqmean_prog(N_blocks, 0.0);
  vector<double> err_prog(N_blocks,0.0);

  for (int i{}; i < N_blocks; i++) {

    for (int j{}; j <= i; j++) { 
      mean_prog[i] += average[j];
      sqmean_prog[i] += square_average[j];
    }
    mean_prog[i]/=(i+1);      
    sqmean_prog[i]/= (i + 1); 
  
      err_prog[i]=error(mean_prog[i], sqmean_prog[i], i); //<x^2> - <x>^2
  }


  // open file
  ofstream out;
  out.open(filename);
  // check 
  if (!out) {
    cerr << "impossibile aprire il file " << filename << endl;
    return;
  }
  //header
  out<<"#BLOCK \t \t AVERAGE \t \t ERROR"<<endl;
  // print results
  for (int i{}; i < N_blocks; i++) {
    out <<i+1<< " \t \t "<< mean_prog[i] << " \t \t " << err_prog[i] << endl;
  }
  out.close();
}
