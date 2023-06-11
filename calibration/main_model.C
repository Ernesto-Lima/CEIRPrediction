#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <math.h>
#include <vector>

#include "main_model.h"

using namespace std;

void RungeKutta(std::vector<double>& QoI,
                double time,
                const double time_step,
                vector<double> Parameters,
                const unsigned int model_n);

void fmodel(std::vector<double> qoi,
            std::vector<double>& auxiliar,
            double rk_time,
            vector<double> Parameters,
            double time,
            const unsigned int model_n);

void main_code(std::vector<double> Parameters,
  	           std::vector<double> time_vec,
  	           std::vector<double> time_hf,
	           std::vector<std::vector<double>>& sol_qoi,
	           std::vector<std::vector<double>>& sol_hf,
	           const unsigned int model_n){
  /////----- Read argumments -----/////
  const double time_step              = 0.01;
  const unsigned int N_time           = time_vec.size()-1;
  unsigned int n_timesteps            = time_vec[N_time]/time_step;
  unsigned int number_equations       = 2;
  double curr_time                    = time_vec[0];
  unsigned int t_step                 = time_vec[0]/time_step;
  /*
  cout << endl;
  cout << "Initial time:     " << time_vec[0] << endl;
  cout << "Final time:       " << time_vec[N_time] << endl;
  cout << "Time to run:      " << time_vec[N_time]-time_vec[0] << endl;
  cout << "Timestep:         " << time_step << endl;
  cout << "Final timestep:   " << n_timesteps << endl;
  cout << "Initial timestep: " << t_step << endl;
  cout << endl;
  */
  /////----- Times to save -----/////
  std::vector<unsigned int> saved_times(N_time+1);
  for(unsigned int i=0; i<=N_time; i++)
    saved_times[i] = time_vec[i]/time_step;
  const unsigned int hf_size = time_hf.size();
  std::vector<unsigned int> hf_times(hf_size);
  for(unsigned int i=0; i<hf_size; i++)
    hf_times[i] = time_hf[i]/time_step;
  /////----- Define initial condition -----/////
  std::vector<double> QoI(number_equations,0);
  QoI[0] = sol_qoi[0][0];
  unsigned int last_p = 0;
  if(model_n >= 10){
    last_p = 8;
  }
  else if(model_n >= 6){
    last_p = 6;
  }
  else if(model_n >= 4){
    last_p = 4;
  }
  else if(model_n >= 2){
    last_p = 3;
  }
  if(model_n >= 2)
    QoI[1] = (1.0-Parameters[last_p]);
  unsigned int pos_sol=1;
  unsigned int pos_hf=0;
  do{
    RungeKutta(QoI,curr_time,time_step,Parameters,model_n);
    t_step++;
    curr_time = t_step*time_step;
    if(t_step==saved_times[pos_sol]){
      sol_qoi[0][pos_sol] = QoI[0];
      pos_sol++;
    }
    if (model_n >= 2){
      if(t_step==hf_times[pos_hf] && t_step<=hf_times[hf_size-1]){
        sol_hf[0][pos_hf] = QoI[1];
        pos_hf++;
      }
    }
  }while(t_step<n_timesteps);
}

void RungeKutta(std::vector<double>& QoI,
                double time,
                const double time_step,
                vector<double> Parameters,
                const unsigned int model_n){

  const unsigned int N = QoI.size();
  std::vector<double> k(4*N,0.0);
  std::vector<double> qoi(N,0.0);
  std::vector<double> auxiliar(N,0.0);
  double rk_time;
  
  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i];}
  rk_time = time;
  fmodel(qoi,auxiliar,rk_time,Parameters,time,model_n);
  for(unsigned int i = 0; i < N; i++){k[i*4+0] = time_step*auxiliar[i];}
  
  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i]+(k[i*4+0]/2.);}
  rk_time = time+0.5*time_step;
  fmodel(qoi,auxiliar,rk_time,Parameters,time,model_n);
  for(unsigned int i = 0; i < N; i++){k[i*4+1] = time_step*auxiliar[i];}
 
  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i]+(k[i*4+1]/2.);}
  rk_time = time+0.5*time_step;
  fmodel(qoi,auxiliar,rk_time,Parameters,time,model_n);
  for(unsigned int i = 0; i < N; i++){k[i*4+2] = time_step*auxiliar[i];}

  for(unsigned int i = 0; i < N; i++){qoi[i] = QoI[i]+k[i*4+2];}
  rk_time = time+time_step;
  fmodel(qoi,auxiliar,rk_time,Parameters,time,model_n);
  for(unsigned int i = 0; i < N; i++){k[i*4+3] = time_step*auxiliar[i];}
  
  for(unsigned int i = 0; i < N; i++){QoI[i] = QoI[i]+(1./6.)*(k[i*4+0]+2*(k[i*4+1]+k[i*4+2])+k[i*4+3]);}

}

void fmodel(std::vector<double> qoi,
            std::vector<double>& auxiliar,
            double rk_time,
            vector<double> Parameters,
            double time,
            const unsigned int model_n){
  double tumor_growth         = Parameters[0];
  if(model_n == 1){
    auxiliar[0] = tumor_growth*qoi[0];
  }
  else if(model_n == 2){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    auxiliar[0] = tumor_growth*qoi[0];
    auxiliar[1] = vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 3){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    auxiliar[0] = tumor_growth*qoi[0]*qoi[1];
    auxiliar[1] = vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 4){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    auxiliar[0] = tumor_growth*qoi[0]
                - immune_response*qoi[0]*qoi[1];
    auxiliar[1] = vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 5){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    auxiliar[0] = tumor_growth*qoi[0]*qoi[1]
                - immune_response*qoi[0]*qoi[1];
    auxiliar[1] = vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 6){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    double gamma                = Parameters[4];
    double r                    = Parameters[5];
    double D                    = 1.0;
    if(time>=7){
      unsigned int n_treat = (time-7)/3;
      for(unsigned int n = 0; n<=n_treat; n++){
        D += gamma * exp(-r*(time-(7+n*3)));
      }
    }
    auxiliar[0] = tumor_growth*qoi[0]
                - D*immune_response*qoi[0]*qoi[1];
    auxiliar[1] = vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 7){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    double gamma                = Parameters[4];
    double r                    = Parameters[5];
    double D                    = 1.0;
    if(time>=7){
      unsigned int n_treat = (time-7)/3;
      for(unsigned int n = 0; n<=n_treat; n++){
        D += gamma * exp(-r*(time-(7+n*3)));
      }
    }
    auxiliar[0] = tumor_growth*qoi[0]*qoi[1]
                - D*immune_response*qoi[0]*qoi[1];
    auxiliar[1] = vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 8){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    double gamma_e              = Parameters[4];
    double r_e                  = Parameters[5];
    double D_e                  = 1.0;
    if(time>=13)
      D_e += gamma_e * exp(-r_e*(time-13));
    if(time>=14)
      D_e += gamma_e * exp(-r_e*(time-14));
    if(time>=15)
      D_e += gamma_e * exp(-r_e*(time-15));
    if(time>=16)
      D_e += gamma_e * exp(-r_e*(time-16));
    if(time>=17)
      D_e += gamma_e * exp(-r_e*(time-17));
    auxiliar[0] = tumor_growth*qoi[0]
                - immune_response*qoi[0]*qoi[1];
    auxiliar[1] = D_e*vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 9){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    double gamma_e              = Parameters[4];
    double r_e                  = Parameters[5];
    double D_e                  = 1.0;
    if(time>=13)
      D_e += gamma_e * exp(-r_e*(time-13));
    if(time>=14)
      D_e += gamma_e * exp(-r_e*(time-14));
    if(time>=15)
      D_e += gamma_e * exp(-r_e*(time-15));
    if(time>=16)
      D_e += gamma_e * exp(-r_e*(time-16));
    if(time>=17)
      D_e += gamma_e * exp(-r_e*(time-17));
    auxiliar[0] = tumor_growth*qoi[0]*qoi[1]
                - immune_response*qoi[0]*qoi[1];
    auxiliar[1] = D_e*vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 10){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    double gamma                = Parameters[4];
    double r                    = Parameters[5];
    double gamma_e              = Parameters[6];
    double r_e                  = Parameters[7];
    double D_e                  = 1.0;
    double D                    = 1.0;
    if(time>=7){
      unsigned int n_treat = (time-7)/3;
      for(unsigned int n = 0; n<=n_treat; n++){
        D += gamma * exp(-r*(time-(7+n*3)));
      }
    }
    if(time>=13)
      D_e += gamma_e * exp(-r_e*(time-13));
    if(time>=14)
      D_e += gamma_e * exp(-r_e*(time-14));
    if(time>=15)
      D_e += gamma_e * exp(-r_e*(time-15));
    if(time>=16)
      D_e += gamma_e * exp(-r_e*(time-16));
    if(time>=17)
      D_e += gamma_e * exp(-r_e*(time-17));
    auxiliar[0] = tumor_growth*qoi[0]
                - D*immune_response*qoi[0]*qoi[1];
    auxiliar[1] = D_e*vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
  else if(model_n == 11){
    double vessel_growth        = Parameters[1];
    double vessel_deterioration = Parameters[2];
    double immune_response      = Parameters[3];
    double gamma                = Parameters[4];
    double r                    = Parameters[5];
    double gamma_e              = Parameters[6];
    double r_e                  = Parameters[7];
    double D_e                  = 1.0;
    double D                    = 1.0;
    if(time>=7){
      unsigned int n_treat = (time-7)/3;
      for(unsigned int n = 0; n<=n_treat; n++){
        D += gamma * exp(-r*(time-(7+n*3)));
      }
    }
    if(time>=13)
      D_e += gamma_e * exp(-r_e*(time-13));
    if(time>=14)
      D_e += gamma_e * exp(-r_e*(time-14));
    if(time>=15)
      D_e += gamma_e * exp(-r_e*(time-15));
    if(time>=16)
      D_e += gamma_e * exp(-r_e*(time-16));
    if(time>=17)
      D_e += gamma_e * exp(-r_e*(time-17));
    auxiliar[0] = tumor_growth*qoi[0]*qoi[1]
                - D*immune_response*qoi[0]*qoi[1];
    auxiliar[1] = D_e*vessel_growth*(1.0-qoi[1])
                - vessel_deterioration*qoi[0]*qoi[1];
  }
}
