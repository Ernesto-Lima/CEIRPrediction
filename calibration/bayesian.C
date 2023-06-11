#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <queso/ConcatenationSubset.h>
#include "GetPot"
#include "main_model.h"

using namespace std;

struct likelihoodRoutine_DataType {
  //---------------------------------------------------------
  //
  // Full_Data is 3D array:
  // 1st dimension: experimental group
  // 2nd dimension: QoIs = [0] time, [1] mean, and [2] standard deviation
  // 3rd dimension: measurement point
  //
  //---------------------------------------------------------
  std::vector< std::vector< std::vector<double> > > TV_Data;
  std::vector< std::vector< std::vector<double> > > HF_Data;
  std::vector<double> Cal_Par;
  unsigned int model_n;
  unsigned int input_n;
  unsigned int position;
  unsigned int N_par;
  unsigned int read_file;
};

double likelihoodRoutine(const QUESO::GslVector& paramValues,
                         const QUESO::GslVector* paramDirection,
                         const void*             functionDataPtr,
                         QUESO::GslVector*       gradVector,
                         QUESO::GslMatrix*       hessianMatrix,
                         QUESO::GslVector*       hessianEffect);
			 
inline void read_full_data(std::vector< std::vector< std::vector<double> > >& v_like_data,
                           std::vector< std::string > name_vector);
                           
int main(int argc, char* argv[]){
  //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  // Starts MPI and checks the input file name
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  MPI_Init(&argc,&argv);
  UQ_FATAL_TEST_MACRO(argc < 2,QUESO::UQ_UNAVAILABLE_RANK,"main()","input file must be argv[1]");
  QUESO::FullEnvironment env(MPI_COMM_WORLD,argv[1],"",NULL);
  //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  // Read model option and file names from file
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  cout << "***************************************************************\n";
  cout << "* Code inputs\n";
  cout << "***************************************************************\n";
  GetPot input_model("model.in");
  const unsigned int model_n = input_model("model_n",1);
  const unsigned int input_n = input_model("input_n",1);
  const unsigned int position = input_model("position",1);
  const unsigned int n_files = input_model("n_files",1);
  const unsigned int read_file = input_model("read_file",0);
  std::vector< std::string > pars_file;
  std::vector< std::string > tv_name;
  for(unsigned int i = 0; i < n_files; i++)
    tv_name.push_back(input_model("all_tv", "nobody",i));
  cout << "Model           = " << model_n << "\n"
       << "Input           = " << input_n << "\n" 
       << "Number of Files = " << n_files << "\n"
       << "TV Files        = " << tv_name[0] << endl;
  for(unsigned int i = 1; i < n_files; i++)
    cout << "                  " << tv_name[i] << endl;
  std::vector< std::string > hf_name;
  for(unsigned int i = 0; i < n_files; i++)
    hf_name.push_back(input_model("all_hf", "nobody",i));
  cout << "HF Files        = " << hf_name[0] << endl;
  for(unsigned int i = 1; i < n_files; i++)
    cout << "                  " << hf_name[i] << endl;
    
  std::vector<double> calibrated_pars;
  if(read_file > 0){
    for(unsigned int i = 0; i < read_file; i++){
        pars_file.push_back(input_model("pars_file", "nobody",i));   
        std::ifstream read;
        read.open(pars_file[i]);
        if(!read.is_open())
          std::cout << "Error opening parameter file.\n";
        std::string line;
        while(std::getline(read, line)){
          double aux;
          std::istringstream(line) >> aux;
          calibrated_pars.push_back(aux);
        }
        read.close();
    }
  }
  for(unsigned int i=0; i<calibrated_pars.size(); i++)
    cout << calibrated_pars[i] << endl;
  //------------------------------------------------------
  //
  //                 CALIBRATION
  //
  //------------------------------------------------------
  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  unsigned int N_par = 1;
  if(model_n == 2 || model_n == 3){
    N_par += 2 + n_files;
  }
  else if(model_n == 4 || model_n == 5){
    N_par += 3 + n_files;
  }
  else if(model_n >= 6 && model_n <= 9){
    N_par += 5 + n_files;
  }
  else if(model_n >= 10){
    N_par += 7 + n_files;
  }
  if(read_file == 1)
    N_par -= calibrated_pars.size()-3;
  if(read_file == 2 && position < 8)
    N_par = 2;
  cout << "Number of pars  = " << N_par << endl; 
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", N_par, NULL);
  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  unsigned int last_p = 0;
  if(read_file == 0){
      paramMins[0] = 0.0;   // k_T
      paramMaxs[0] = 0.4; //0.4
      if(model_n >= 2){
        paramMins[1] = 0.0; // k_V
        paramMaxs[1] = 0.5; //0.5
        paramMins[2] = 0.0; // \mu_V
        paramMaxs[2] = 0.5; //0.5
        last_p = 3;
      } 
      if(model_n >= 4){
        paramMins[3] = 0.0; // mu_T
        paramMaxs[3] = 0.5;
        last_p = 4;
      }
      if(model_n >= 6){
        paramMins[4] = 0.0;  // \gamma
        paramMaxs[4] = 3.5; //50
        paramMins[5] = 0.0;  // r
        paramMaxs[5] = 0.2; //10
        last_p = 6;
      }
      if(model_n >= 10){
        paramMins[6]  = 0.0;  // \gamma_e
        paramMaxs[6]  = 3.5;
        paramMins[7] = 0.0;  // r_e
        paramMaxs[7] = 0.2;
        last_p = 8;
      }
      if(model_n >= 2){
        for(unsigned int p = 0; p<n_files; p++){
          paramMins[last_p+p] = 0.0; // H_f
          paramMaxs[last_p+p] = 1.0;
        }
      }
  }
  else if(read_file == 1){
      if(model_n >= 10){
        paramMins[0]  = 0.0;  // \gamma_e
        paramMaxs[0]  = 3.5;
        paramMins[1] = 0.0;  // r_e
        paramMaxs[1] = 0.2;
        last_p = 2;
      }
      if(model_n >= 2){
        for(unsigned int p = 0; p<n_files; p++){
          paramMins[last_p+p] = 0.0; // H_f
          paramMaxs[last_p+p] = 1.0;
        }
      }
  }
  else if(read_file == 2){
    paramMins[0] = 0.0; // H_f
    paramMaxs[0] = 1.0;
    if(position == 0){
      paramMins[1] = 0.0;   // k_T
      paramMaxs[1] = 0.4; //0.4
    }
    else if(position == 1){
      paramMins[1] = 0.0; // k_V
      paramMaxs[1] = 0.5; //0.5
    }
    else if(position == 2){
      paramMins[1] = 0.0; // \mu_V
      paramMaxs[1] = 0.5; //0.5
    }
    else if(position == 3){
      paramMins[1] = 0.0; // mu_T
      paramMaxs[1] = 0.5;
    }
    else if(position == 4){
      paramMins[1] = 0.0;  // \gamma
      paramMaxs[1] = 3.5; //50
    }
    else if(position == 5){
      paramMins[1] = 0.0;  // r
      paramMaxs[1] = 0.2; //10
    }
    else if(position == 6){
      paramMins[1]  = 0.0;  // \gamma_e
      paramMaxs[1]  = 3.5;
    }
    else if(position == 7){
      paramMins[1] = 0.0;  // r_e
      paramMaxs[1] = 0.2;
    }
  }
  for(unsigned int i=0; i<N_par; i++)
    cout << "\t\t\t\t\t Min and max = " << paramMins[i] << ", " << paramMaxs[i] << endl;
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",
                paramSpace,
                paramMins,
                paramMaxs);
  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function
  //------------------------------------------------------
  cout << "***************************************************************\n";
  cout << "* Measured times from each data point\n";
  cout << "***************************************************************\n";
  likelihoodRoutine_DataType likelihoodRoutine_Data;
  read_full_data(likelihoodRoutine_Data.TV_Data,tv_name);
  read_full_data(likelihoodRoutine_Data.HF_Data,hf_name);
  likelihoodRoutine_Data.N_par = N_par;
  cout << "Calibrating the model using " << likelihoodRoutine_Data.TV_Data.size() << " data sets\n";
  for (unsigned int m=0; m<likelihoodRoutine_Data.TV_Data.size(); m++){
    for (unsigned int i=0; i<likelihoodRoutine_Data.TV_Data[m][0].size(); i++)
      cout << likelihoodRoutine_Data.TV_Data[m][0][i] << " ";
    cout << "<-- Measured times\n";
    for (unsigned int i=0; i<likelihoodRoutine_Data.HF_Data[m][0].size(); i++)
      cout << likelihoodRoutine_Data.HF_Data[m][0][i] << " ";
    cout << "<-- Measured times\n\n";
  }
  likelihoodRoutine_Data.model_n = model_n;
  likelihoodRoutine_Data.input_n = input_n;
  likelihoodRoutine_Data.position = position;
  likelihoodRoutine_Data.read_file = read_file;
  likelihoodRoutine_Data.Cal_Par = calibrated_pars;
  for(unsigned int i=0; i<likelihoodRoutine_Data.Cal_Par.size(); i++)
    cout << likelihoodRoutine_Data.Cal_Par[i] << endl;
  //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  // Generate generic scalar function
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("cal_like_",
	                      paramDomain,
	                      likelihoodRoutine,
	                      (void *) &likelihoodRoutine_Data,
	                      true);
  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix> priorRv("cal_prior_", paramDomain);
  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> postRv("cal_post_", paramSpace);
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("",
	   NULL,
	   priorRv,
	   likelihoodFunctionObj,
	   postRv);
  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV
  //------------------------------------------------------
  cout << "***************************************************************\n";
  cout << "* Starting multi-level sampling\n";
  cout << "***************************************************************\n";
  ip.solveWithBayesMLSampling();
  MPI_Finalize();
  return 0;
}
//------------------------------------------------------
//
// Inverse Problem Routine
//
//------------------------------------------------------
double likelihoodRoutine(const QUESO::GslVector& paramValues,
                         const QUESO::GslVector* paramDirection,
                         const void*             functionDataPtr,
                         QUESO::GslVector*       gradVector,
                         QUESO::GslMatrix*       hessianMatrix,
                         QUESO::GslVector*       hessianEffect){
  //cout << "\t\t\t\t\t Okay!"<< endl;
  //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  // Model Parameters
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //cout << "\t\t\t\t\t Read pars!"<< endl;
  const vector< vector< vector<double> > >
    & TV_Data = ((likelihoodRoutine_DataType *) functionDataPtr)->TV_Data;
  //cout << "\t\t\t\t\t TV Data!"<< endl;
  const vector< vector< vector<double> > >
    & HF_Data = ((likelihoodRoutine_DataType *) functionDataPtr)->HF_Data;
  //cout << "\t\t\t\t\t HF Data!"<< endl;
  const vector<double>
    & Cal_Par = ((likelihoodRoutine_DataType *) functionDataPtr)->Cal_Par;
  const unsigned int model_n = ((likelihoodRoutine_DataType *) functionDataPtr)->model_n;
  const unsigned int input_n = ((likelihoodRoutine_DataType *) functionDataPtr)->input_n;
  const unsigned int position = ((likelihoodRoutine_DataType *) functionDataPtr)->position;
  const unsigned int read_file = ((likelihoodRoutine_DataType *) functionDataPtr)->read_file;
  //cout << "\t\t\t\t\t Model!"<< endl;
  unsigned int N_par = ((likelihoodRoutine_DataType *) functionDataPtr)->N_par;
  if (read_file == 1)
    N_par += Cal_Par.size()-3;
  else if(read_file == 2)
    N_par = 9;
  //cout << N_par << "\t\t\t\t\t Par!"<< endl;
  const QUESO::BaseEnvironment& env = paramValues.env();
  //cout << "\t\t\t\t\t param values!"<< endl;
  std::vector<double> Parameters(N_par,0.0);
  //cout << "\t\t\t\t\t assigned par!"<< endl;
  unsigned int n_qois = 1;
  double result_t = 0.;
  double result_v = 0.;
  const unsigned int data_size = TV_Data.size();
  //cout << "\t\t\t\t\t Data size!"<< endl;
  for (unsigned int s = 0; s < data_size; s++){
    if(read_file == 1){
        for (int i=0; i<int(N_par); i++){
            if(i<int(Cal_Par.size())-3)
                Parameters[i] = Cal_Par[i];
            else
                Parameters[i] = paramValues[i-int(Cal_Par.size())+3];
        }
    }
    else if(read_file == 2){
        for (int i=0; i<int(N_par); i++){
            if(i<6)
                Parameters[i] = Cal_Par[i];
            if(i<8)
                Parameters[i] = Cal_Par[i+3];
            else
                Parameters[i] = paramValues[i-8];
        }
        if (position < 8)
          Parameters[position] = paramValues[1];
    }
    else{
        for (int i=0; i<int(N_par); i++)
            Parameters[i] = paramValues[i];
    }
    if(model_n >= 6 && input_n == 1 && s == 0)
        Parameters[4] = 0.0;
    if(model_n >= 10 && input_n == 4 && s == 0)
        Parameters[6] = 0.0;
    //cout << "\t\t\t\t\t Looping starting!"<< endl;
    unsigned int tv_size = TV_Data[s][0].size();
    unsigned int hf_size = HF_Data[s][0].size();
    //cout << "\t\t\t\t\t data!"<< endl;
    std::vector<std::vector<double>> sol_tv(n_qois, std::vector<double> (tv_size,0));
    std::vector<std::vector<double>> sol_hf(n_qois, std::vector<double> (hf_size,0));
    //cout << "\t\t\t\t\t vector!"<< endl;
    //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    // Run model
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    sol_tv[0][0] = TV_Data[s][1][0];
    //cout << "\t\t\t\t\t ic! " << N_par-data_size-1 << ", " << N_par-data_size-1+s << endl;
    //cout << "Npar = " << N_par << ", datasize = " << data_size << endl;
    if(model_n>=2){
        if(read_file == 1)
            Parameters[N_par-data_size] = paramValues[(N_par-(Cal_Par.size()-3))-data_size+s];
        else if(read_file == 2)
            Parameters[N_par-data_size] = paramValues[0];
        else
            Parameters[N_par-data_size] = paramValues[N_par-data_size+s];
    }
    //cout << "Let's solve it" << endl;
    main_code(Parameters,TV_Data[s][0],HF_Data[s][0],sol_tv,sol_hf,model_n);
    //cout << "Solved!" << endl;
    if(env.subRank() == 0){
      double total_tv = 0.0;
      double total_hf = 0.0;
      double res_tv   = 0.0;
      double res_hf   = 0.0;
      for(unsigned int i=1; i<tv_size; i++){
        if(model_n >= 2){
          for(unsigned int j=0; j<hf_size; j++){
            if(TV_Data[s][0][i] == HF_Data[s][0][j]){
              double sol_model = 1.0-sol_hf[0][j];
              double sol_data  = HF_Data[s][1][j];
              double std_dev   = HF_Data[s][2][j];
              double ratio_v   = (sol_model - sol_data)/std_dev;
              total_hf += sol_data;
              res_hf += ratio_v*ratio_v
                      + log(std_dev)
                      + 0.5*log(M_PI)
                      + 0.5*log(2.0);
            }
          }
        }
        double sol_model = sol_tv[0][i];
        double sol_data  = TV_Data[s][1][i];
        double std_dev   = TV_Data[s][2][i];
        double ratio_t   = (sol_model - sol_data)/std_dev;        
        total_tv += sol_data;
        res_tv += ratio_t*ratio_t
                + log(std_dev)
                + 0.5*log(M_PI)
                + 0.5*log(2.0);
      }
      if(model_n >= 2)
        result_v += res_hf;///total_hf;
      result_t += res_tv;///total_tv;
    }
  }
  //cout << result_t << " " << result_v << endl;
  //getchar();
  return -100.0*(1.0*result_t+1.0*result_v);
}

inline void read_full_data(std::vector< std::vector< std::vector<double> > >& v_like_data,
                           std::vector< std::string > name_vector){
  const unsigned int size = name_vector.size();
  for(unsigned int f = 0; f<size; f++){
    std::vector< std::vector<double> > data_vec;
    std::vector<double> aux_vec_0;
    std::vector<double> aux_vec_1;
    std::vector<double> aux_vec_2;
    std::ifstream read;
    read.open(name_vector[f]);
    if(!read.is_open())
      std::cout << "Error opening data file " << name_vector[f] << ".\n";
    std::string line;
    while(std::getline(read, line)){
      double aux[3];
      std::istringstream(line) >> aux[0] >> aux[1] >> aux[2];
      aux_vec_0.push_back(aux[0]);
      aux_vec_1.push_back(aux[1]);
      aux_vec_2.push_back(aux[2]);
    }
    data_vec.push_back(aux_vec_0);
    data_vec.push_back(aux_vec_1);
    data_vec.push_back(aux_vec_2);
    v_like_data.push_back(data_vec);
    read.close();
  }
}
