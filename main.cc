#include<algorithm>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<iomanip>
#include<map>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>

#include<unistd.h>
#include<getopt.h> 

#include"matrix/matrix.h"
#include"matrix/myVector.h"
#include"distribution/distribution.h"

#define EPSILON_H 1e-14
#define EPSILON_G 1e-10
#define EPSILON_Q 1e-11
#define EPSILON_V 1e-10
#define EPSILON_Y 1e-14
#define EPSILON_B_SUM 1e-16
#define EPSILON_B 1e-14
#define EPSILON_BINOM 1e-13

#define MAX_ITER_Q 10000
#define MAX_ITER 10000

#define ARG_SAMPLING_MEAN  1
#define ARG_SAMPLING_CV    2
#define ARG_SRVMEAN        3
#define ARG_SRVCV          4
#define ARG_LMD_BG         5
#define ARG_SRVMEAN_BG     6
#define ARG_SRVCV_BG       7
#define ARG_MU             8
#define ARG_OUTPUT         9

void read_arg(int argc, char* argv[],  PhDist*& SamplingDist, Distribution*& srvDist, double& lmd_bg, Distribution*& srvDist_bg, double& mu, std::string& out_file);

// Compute the invariant probability vector of a transition rate matrix by iteration
RowVector computeIPV(const Matrix& mat, double tolerance, int max_iter);

int main(int argc, char* argv[]){
  // Set output precision
  std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(14);

  // Output filename
  std::string out_file;

  // Arrival rate of the background stream
  double lmd_bg = 0.0;

  // Service time distribution of the background stream
  Distribution* srvDist_bg_ptr;

  // Service time distribution of the tagged stream
  Distribution* srvDist_ptr;
	
  // Inter-arrival time distribution (tagged stream)
  PhDist* arvDist_ptr;

  // Service rate
  double mu = 0.0;

  read_arg(argc, argv, arvDist_ptr, srvDist_ptr, lmd_bg, srvDist_bg_ptr, mu, out_file);

  const PhDist& arvDist = *arvDist_ptr;
  const Distribution& srvDist = *srvDist_ptr;
  const Distribution& srvDist_bg = *srvDist_bg_ptr;

  double EG = arvDist._mean;
  double CvG = arvDist._coef_v;
  double lmd = 1.0/EG;

  RowVector gamma = arvDist._alpha;
  Matrix Gamma = arvDist._T;
  int M = Gamma.size();

  double EH = srvDist._mean;
  double CvH = srvDist._coef_v;

  double EH_bg = srvDist_bg._mean;
  double CvH_bg = srvDist_bg._coef_v;

  double rho_bg = lmd_bg*srvDist_bg._mean/mu;
  double rho = lmd*srvDist._mean/mu;

  Matrix I(M,M);
  I.setIdentity();

  ColVector e(M);
  e.setAllOne();

  std::cout << "----------------------------------------------\n";
  std::cout << "mu=" << mu << std::endl;
  std::cout << "lambda_bg=" << lmd_bg << std::endl;
  std::cout << "Service time distribution (background stream):" << std::endl;
  srvDist_bg.print();
  std::cout << std::endl;

  std::cout << "Inter-Sampling time distribution:" << std::endl;
  arvDist.print();
  std::cout << std::endl;

  std::cout << "lambda=" << lmd << std::endl;
  std::cout << "Service time distribution (tagged stream):" << std::endl;
  srvDist.print();
  std::cout << "----------------------------------------------\n";

  double theta = Gamma.getMaxAbsDiag();
  double zeta = theta + lmd_bg;

  Matrix C = Gamma - I*lmd_bg;

  std::cout << "Computing c" << std::endl;
  ColVector c = arvDist.unifDensityVector(EPSILON_G, theta);

  // Compute hat{h} and hat{h}_bg
  std::cout << "Computing hat_h" << std::endl;
  RowVector hat_h = srvDist.PoissonVec(EPSILON_V*1e-4, zeta/mu);

  std::cout << "Computing hat_h_bg" << std::endl;
  RowVector hat_h_bg = srvDist_bg.PoissonVec(EPSILON_V*1e-4, zeta/mu);

  int hat_h_size = hat_h.size();
  int hat_h_bg_size = hat_h_bg.size();

  int max_h_size = std::max(hat_h_size, hat_h_bg_size);

  // Compute hatD;
  std::vector<Matrix> hatD; 
  for(int m=0; m < max_h_size; m++){
    // Initialize by a zero matrix
    Matrix hatD_m(M,M);

    if(hat_h_size > m){
      hatD_m = hatD_m + hat_h[m]*(-1.0)*Gamma*e*gamma;
    }

    if(hat_h_bg_size > m){
      hatD_m = hatD_m + lmd_bg*hat_h_bg[m]*I;
    }

    hatD.push_back(hatD_m);
  }

  Matrix Q = C;

  std::cout << "#Computing Q" << std::endl;
  for(int i=0; i < MAX_ITER; i++){
    Matrix prevQ = Q;
    Q = C;

    Matrix P = I;

    for(int k=0; k < hatD.size(); k++){
      Q = Q + hatD[k]*P;
      P = P * (I + (1.0/zeta)*prevQ);
    }

    double residual = (Q*e).getMaxAbs();

    std::cout << "max|Qe|=" << residual << std::endl;
    if(residual < EPSILON_Q){
      break;
    }

    if(i == MAX_ITER -1){
       std::cerr << "Computation of Q not convergent" << std::endl;
       return 1;
    }
  }

  ColVector Qe = Q*e;
  for(int i=0; i < Q.size(); i++){
    Q[i][i] -= Qe[i];
  }

  std::cout << "#Computing kappa" << std::endl;
  RowVector kappa = computeIPV(Q, EPSILON_V, MAX_ITER);

  std::vector<Matrix> R;

  std::cout << "#Computing R" << std::endl;
  for(int m=0; m+1 < hatD.size(); m++){
    Matrix R_m(M,M);
    Matrix P = I;

    for(int k=0; m+k+1 < hatD.size(); k++){
      R_m = R_m + (1.0/zeta)*hatD[m+k+1]*P;
      P = P*(I + (1.0/zeta)*Q);
    }

    R.push_back(R_m);
  }

  Matrix R_sum(M,M);
  for(int i=0; i < R.size(); i++){
    R_sum = R_sum + R[i];
  }

  std::vector<RowVector> hat_v;
  std::cout << "#Computing hat{v}" << std::endl;

  {
    double v_sum = 0.0;
    Matrix invR0 = (I - R[0]).getInverse();
     
    RowVector hat_v_0 = (kappa*invR0)*(1.0 - rho - rho_bg);
    hat_v.push_back(hat_v_0);

    v_sum += hat_v_0*e;

    for(int m=1; m < MAX_ITER; m++){
      RowVector hat_v_m(M);

      for(int k=0; k <= m-1; k++){
        if(m-k < R.size()){
          hat_v_m = hat_v_m + hat_v[k]*R[m-k];
        }
      }

      hat_v_m = hat_v_m*invR0;
      hat_v.push_back(hat_v_m);

      v_sum += hat_v_m*e;

      if(v_sum > 1.0 - EPSILON_V){
         break;
      }

      if(m == MAX_ITER -1){
         std::cerr << "Computation of v not convergent" << std::endl;
         return 1;
      }
    }
  }

  RowVector hat_v_sum(M);
  for(int i=0; i < hat_v.size(); i++){
    hat_v_sum = hat_v_sum + hat_v[i];
  }

  std::vector<double> d;
  std::cout << "#Computing d" << std::endl;

  double ED = 0.0;

  {
    double d_sum = 0.0;
    for(int m=0; m+1 < hat_v.size() + hat_h.size(); m++){
      double d_m = 0.0;

      for(int k=0; k <= m; k++){
        if( ( m-k < hat_h.size() ) && k < hat_v.size()){
          d_m += (1.0/lmd)*(hat_v[k]*((-1.0)*Gamma*e))*hat_h[m-k];
        }
      }

      d_sum += d_m;
      ED += (1.0/zeta)*m*d_m;

      d.push_back(d_m);

      if(d_sum > 1.0 - EPSILON_V){
        break;
      }
    }
  }

  std::vector< std::vector<double> > y;
  
  std::vector<double> y0;
  std::vector<int> y_s;
  y0.push_back(1.0);
  y.push_back(y0);

  // We only need to compute b[i] up to i = c.size()
  std::vector<double> b(c.size()+1, 0.0);
  double b_sum = 0.0;
  std::cout << "#Computing b" << std::endl;

  for(int k=0; k < MAX_ITER; k++){
    std::vector<double> y_k1;
    double y_sum = 0.0;

    for(int m=0; m < hat_h_bg.size() + y[k].size(); m++){
      double y_val = 0.0;
      for(int i=0; i <= m; i++){

        if( (y[k].size() > i) && (hat_h_bg.size() > m-i) ){
          y_val += y[k][i] * hat_h_bg[m-i];
        }
      }

      y_k1.push_back(y_val);
      y_sum += y_val;

      if(y_sum > 1.0 - EPSILON_Y){
        break;
      }
    }

    y.push_back(y_k1);

    int binom_start = 0;
    int binom_end = 0;
    std::vector<double> binom = BinomialFunc(k, theta/(theta+lmd_bg), EPSILON_BINOM, binom_start, binom_end);

    double b_sum_delta = 0.0;

    for(int i=0; i <= std::min(k, int(b.size())-1); i++){

      if(k < y[k-i+1].size() && i >= binom_start && i <= binom_end){
        double delta = (1.0/(k-i+1.0)) * binom[i] * y[k-i+1][k];
        b[i] += delta;
        b_sum_delta += delta;
      }
    }

    b_sum += b_sum_delta;
    std::cout << "k = " << k << ", b_sum = " << b_sum << std::endl;

    if(b_sum_delta < EPSILON_B_SUM){
      break;
    }

    if(k == MAX_ITER -1){
       std::cerr << "Computation of b not convergent" << std::endl;
       return 1;
    }
    
  }

  for(int i=0; i < b.size(); i++){
    b[i] = b[i]*lmd_bg /(theta + lmd_bg);
  }

  b[1] += theta/(theta + lmd_bg);

  std::vector<double> q;
  std::vector< std::vector<double> > b_conv;

  for(int l=0; l < d.size(); l++){
  }
  std::cout << "Computing b_l" << std::endl;

  std::vector<double> b_0(b.size(), 0.0); 
  b_0[0] = 1.0;
  b_conv.push_back(b_0);

  for(int l=1; l < MAX_ITER; l++){
    std::vector<double> b_l(b.size(), 0.0); 
    double b_l_sum = 0.0;

    for(int m=0; m < b.size(); m++){
      for(int i=0; i <= m; i++){
        b_l[m] += b_conv[l-1][i] * b[m-i];
        b_l_sum += b_l[m];
      }
    }

    b_conv.push_back(b_l);

    if(b_l_sum < EPSILON_B){
      break;
    }

    if(l == MAX_ITER -1){
       std::cerr << "Computation of b_l not convergent" << std::endl;
       return 1;
    }
  }

  std::vector<double> d_sum(d.size(), 0.0);
  d_sum[0] = d[0];
  for(int m=1; m < d_sum.size(); m++){
    d_sum[m] = d_sum[m-1] + d[m];
  }

  double q_0 = ED/theta - (1.0-rho_bg)/(theta*theta);
  for(int l=0; l < b_conv.size(); l++){
    double delta_q_0 = 0.0;

    if(l < d_sum.size()){
      delta_q_0 = d_sum[l];

    } else{
      delta_q_0 = 1.0;
    }

    delta_q_0 = (1.0/(theta*zeta))*b_conv[l][0]*delta_q_0;
    q_0 = q_0 + delta_q_0;
  }

  q.push_back(q_0);

  double PK_formula = 0.5*rho_bg*(1.0+CvH_bg*CvH_bg)/(1.0-rho_bg)/theta/mu;


  std::cout << "q[" << 0 << "] = " << q[0] << std::endl;

  for(int k=1; k < c.size()+1; k++){
    double q_k = q[k-1] - (1.0-rho_bg)/(theta*theta);

    for(int l=0; l < b_conv.size(); l++){
      double delta_q_k = 0.0;

      if(l < d_sum.size()){
        delta_q_k = d_sum[l];

      }else{
        delta_q_k = 1.0;
      }

      delta_q_k = (1.0/(theta*zeta))*b_conv[l][k]*delta_q_k;
      q_k = q_k + delta_q_k;
    }

    q.push_back(q_k);
    std::cout << "q[" << k << "] = " << q[k] << std::endl;
  }

  std::cout << "PK-formula: " << PK_formula  << std::endl;

  double EA = EH/mu;

  for(int k=0; k < c.size(); k++){
    EA += (c[k]/(theta*EG))*( (k+1.0)*q[k+1] + (k+1.0)*(k+2.0)/(2.0*theta*theta) );
  }

  std::cout << std::setprecision(8);
  std::cout << "E[A] = " << EA << std::endl;

  std::ofstream res(out_file, std::ios::app);
  res << std::setiosflags(std::ios::fixed) 
      << std::setprecision(8);

  res << lmd << " " << CvG << " " << EH << " " << CvH << " " << lmd_bg << " " << EH_bg << " " << CvH_bg << " " << EA  << " " << mu << std::endl;

  return 0;
}

RowVector computeIPV(const Matrix& gen, double tolerance, int max_iter){
	// Uniformization
	double theta = gen.getMaxAbsDiag()*1.1;
  int M = gen.size();

  if(M == 1){
    RowVector res(1);
    res[0] = 1.0;
    return res;
  }

  Matrix I(M, M);
  I.setIdentity();

  ColVector e(M);
  e.setAllOne();

  Matrix mat = I + gen*(1.0/theta);

  Matrix pow_mat = mat;
  bool isConvergent = false;

  for(int i = 0; i < max_iter; i++){
    pow_mat = pow_mat * pow_mat;

    if( (pow_mat*gen).getMaxAbs() < tolerance){
      isConvergent = true;
      break;
    } 
  }
  if(isConvergent){
    return pow_mat[0]*(1.0/(pow_mat[0]*e));
  }else{
    std::cerr << "ERROR (getIPV): Iteration not convergent" << std::endl;
    exit(1);
  }
}


void read_arg(int argc, char* argv[],  PhDist*& SamplingDist, Distribution*& srvDist, double& lmd_bg, Distribution*& srvDist_bg, double &mu, std::string& out_file){

  struct option long_options[] = {
    {"EG",     required_argument, NULL, ARG_SAMPLING_MEAN},
    {"CvG",    required_argument, NULL, ARG_SAMPLING_CV},
    {"EH",     required_argument, NULL, ARG_SRVMEAN},
    {"CvH",    required_argument, NULL, ARG_SRVCV},
    {"lmd_bg",  required_argument, NULL, ARG_LMD_BG},
    {"EH_bg",  required_argument, NULL, ARG_SRVMEAN_BG},
    {"CvH_bg", required_argument, NULL, ARG_SRVCV_BG},
    {"mu", required_argument, NULL, ARG_MU},
    {"output", required_argument, NULL, ARG_OUTPUT},
    {0, 0, 0, 0}
};

  std::stringstream opt_str;
  opt_str << ARG_SAMPLING_MEAN  << ":"
          << ARG_SAMPLING_CV    << ":"
	  << ARG_SRVMEAN        << ":"
	  << ARG_SRVCV          << ":"
	  << ARG_LMD_BG         << ":"
	  << ARG_SRVMEAN_BG     << ":"
	  << ARG_SRVCV_BG       << ":"
	  << ARG_MU             << ":"
	  << ARG_OUTPUT         << ":";

  double EG = -1;
  double CvG = -1;
  double EH = -1;
  double CvH = -1;

  double EH_bg = -1;
  double CvH_bg = -1;

  int count = 0;

  while(1){
    int option_index = 0;
    int i = getopt_long(argc, argv, opt_str.str().c_str(), long_options, &option_index);

    if(i == -1 && count == 9){break;}

    switch(i){
	case ARG_SAMPLING_MEAN:
		{
			std::stringstream buff(optarg);
			buff >> EG;
        		count += 1;
			break;
		}
	case ARG_SAMPLING_CV: 
		{
			std::stringstream buff(optarg);
			buff >> CvG;
          		count += 1;
			break;
		}
	case ARG_SRVMEAN:
		{
			std::stringstream buff(optarg);
			buff >> EH;
        		count += 1;
			break;
		}
	case ARG_SRVCV: 
		{
			std::stringstream buff(optarg);
			buff >> CvH;
          		count += 1;
			break;
		}
	case ARG_LMD_BG: 
		{
			std::stringstream buff(optarg);
			buff >> lmd_bg;
          		count += 1;
			break;
		}
	case ARG_SRVMEAN_BG:
		{
			std::stringstream buff(optarg);
			buff >> EH_bg;
        		count += 1;
			break;
		}
	case ARG_SRVCV_BG: 
		{
			std::stringstream buff(optarg);
			buff >> CvH_bg;
          		count += 1;
			break;
		}
	case ARG_MU: 
		{
			std::stringstream buff(optarg);
			buff >> mu;
        		count += 1;
			break;
		}
	case ARG_OUTPUT:
		{
			std::stringstream buff(optarg);
			buff >> out_file;
        		count += 1;
			break;
		}
	default:
		{
			std::cout 
			<< "\n"
			<< "Usage:\n"
			<< "    --EG:\t\t" 
			<< "Mean inter-sampling time\n"
			<< "    --CvG:\t\t" 
			<< "CV of inter-sampling time\n"
			<< "    --EH:\t\t" 
			<< "Mean service time of tagged stream\n"
			<< "    --CvH:\t\t" 
			<< "Cv of service time of tagged stream\n"
			<< "    --lmd_bg:\t\t" 
			<< "Arrival rate of background stream\n"
			<< "    --EH_bg:\t\t" 
			<< "Mean service time of background stream\n"
			<< "    --CvH_bg:\t\t" 
			<< "Cv of service time of background stream\n"
			<< "    --mu:\t\t" 
			<< "Service rate\n"
			<< std::endl;
			
			exit(0);
		}
	}
　}

  if(CvG > 0 && CvG < 1){
    SamplingDist = new MixErDist(EG, CvG);

  } else if(CvG == 1.0){
    SamplingDist = new ExpDist(EG);

  } else if(CvG > 1){
    SamplingDist = new BalHyperDist(EG, CvG);

  } else{
    std::cerr << "ERROR (read_arg): Invalid value for the coefficient of variation of SamplingDist" << std::endl;
    exit(1);
  }

  if(CvH == 0){
    srvDist = new ConstDist(EH);

  } else if(CvH > 0 && CvH < 1){
    srvDist = new MixErDist(EH, CvH);

  } else if(CvH == 1.0){
    srvDist = new ExpDist(EH);

  } else if(CvH > 1){
    srvDist = new BalHyperDist(EH, CvH);

  } else{
    std::cerr << "ERROR (read_arg): Invalid value for the coefficient of variation of srvDist" << std::endl;
    exit(1);
  }

  if(CvH_bg == 0){
    srvDist_bg = new ConstDist(EH_bg);

  } else if(CvH_bg > 0 && CvH_bg < 1){
    srvDist_bg = new MixErDist(EH_bg, CvH_bg);

  } else if(CvH_bg == 1.0){
    srvDist_bg = new ExpDist(EH_bg);

  } else if(CvH_bg > 1){
    srvDist_bg = new BalHyperDist(EH_bg, CvH_bg);

  } else{
    std::cerr << "ERROR (read_arg): Invalid value for the coefficient of variation of srvDist" << std::endl;
    exit(1);
  }

}

