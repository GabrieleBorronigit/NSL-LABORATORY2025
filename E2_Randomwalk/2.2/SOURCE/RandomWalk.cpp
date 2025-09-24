#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "RandomWalk.h"
using namespace std;




void RandomWalk:: History_ToFile(string filename){
    //apro il file
    ofstream out;
    out.open(filename);
    //controllo
    if(!out)
    {
        cerr<<"impossibile creare il file "<<filename<<endl;
        return;
    }

//stampo tutta la storia

for(int i{}; i<m_History.size();i++){
    for(int j{};j<m_D;j++)
    {
        out<<m_History[i][j]<<"\t";

    }
    out<<endl;
}
out.close();


}





double RandomWalk::Distance(){
    //distanza euclidea m_D-dimensionale
double d{};
for(int i{};i<m_D;i++)
{
    d+=m_r[i]*m_r[i];
}
return d;

}


void RandomWalk::LatticeWalk(bool save_history){

 

    for(int i{};i<m_N_steps;i++)
    {
        //choose a random discrete direction
        int direction=static_cast<int>(m_gen.Rannyu(0.0,m_D));

    //coin flip for the next step
   int  s= (m_gen.Rannyu()>0.5) ? -1 : +1;
    m_r[direction]+= s*m_a;

    m_amsd_lattice.push_back(Distance());
    if(save_history==true)
    Save_History();
    }

}





//samples a random point on the unitary sphere
vector<double> RandomWalk::Uniform_continuum_direction(){
    //uniform phi in [0,2pi]
    double phi{m_gen.Rannyu(0.0, 2*M_PI)};


    //sample theta with the inverse cumulative function

    double theta{acos(1-2.0*m_gen.Rannyu())};

    vector<double> unit_dir; //point on the unitary sphere
    unit_dir.push_back(sin(theta)*cos(phi)); //x
    unit_dir.push_back(sin(theta)*sin(phi)); //y
    unit_dir.push_back(cos(theta));          //z

return unit_dir;


}

void RandomWalk:: ContinuumWalk(bool save_history){


//loop over steps
for(int i{};i<m_N_steps;i++)
{
    //sample direction
vector<double> direction=Uniform_continuum_direction();

//coin flip
int s= (m_gen.Rannyu()>0.5) ? -1 : +1;

//update each position vector's coordinate
for(int j{};j<m_D;j++){
    m_r[j]+=s*m_a*direction[j];
}

m_amsd_lattice.push_back(Distance());

if(save_history==true)
Save_History();
}



}




double error(double av, double av2, int n) {  // mean dev std

  if (n == 0)
    return 0;  // set first block error to zero
   

  return sqrt((av2 - pow(av, 2)) / n);
}



//computes statistical averages and errors using data blocking, outputs on file
void ComputeMeanErrors_DB(const vector<double> &data, int N_blocks, ofstream&out, int step) {

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

  // print results
  
    out <<step+1<< " \t \t "<< sqrt(mean_prog.back()) << " \t \t " << 0.5*err_prog.back()/sqrt(mean_prog.back()) << endl;
  

}
