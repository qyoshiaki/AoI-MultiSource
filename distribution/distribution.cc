#include<iostream>
#include"distribution.h"

Distribution::Distribution(std::string tp, double m, double cv) 
					: _type(tp), _mean(m), _coef_v(cv)
{
	if(m <= 0){
		std::cerr << "ERROR (Distribution::Distribution): Invalid value for mean" 
				<< std::endl;
		exit(1);
	}
	if(cv < 0){
		std::cerr << "ERROR (Distribution::Distribution): Invalid value for coefficient of variation" 
				<< std::endl;
		exit(1);
	}
}

Distribution::Distribution(const char* tp, double m, double cv)
					: _type(tp), _mean(m), _coef_v(cv){

	if(m <= 0){
		std::cerr << "ERROR (Distribution::Distribution): Invalid value for mean" 
				<< std::endl;
		exit(1);
	}
	if(cv < 0){
		std::cerr << "ERROR (Distribution::Distribution): Invalid value for coefficient of variation" 
				<< std::endl;
		exit(1);
	}
}

void Distribution::print() const{
	std::cout << "#type: " << _type << std::endl;
	std::cout << "#mean: " << _mean << std::endl;
	std::cout << "#coef_v: " << _coef_v << std::endl;
}

void ErDist::print() const{
	std::cout << "#type: " << _type << "_" << _stage << std::endl;
	std::cout << "#mean: " << _mean << std::endl;
	std::cout << "#coef_v: " << _coef_v << std::endl;
}

PhDist::PhDist(RowVector vec, Matrix mat) : _alpha(vec), _T(mat)
{
	_type = std::string("PH_General");

	Matrix invT = (-1)*_T.getInverse();
	ColVector e(_T.size());
	e.setAllOne();

	_mean = _alpha * invT * e;
	_coef_v = sqrt( 2.0*(_alpha*invT*invT*e)/(_mean*_mean) - 1.0 );
}

ErDist::ErDist(double m, int k) : PhDist("Er", m, sqrt(1.0/(double)k)), _stage(k){
	if(k <= 0){
		std::cerr << "ERROR (ErDist::ErDist): Invalid value for number of stages" 
				<< std::endl;
		exit(1);
	}

	double mu = (double)k/m;

	// Initialize with a zero vector
	_alpha = RowVector(k);

	_alpha[0] = 1.0;
	
	// Initialize with a zero matrix
	_T = Matrix(k,k);

	for(int i=0; i < k-1; i++){
		_T[i][i] = (-1.0)*mu;
		_T[i][i+1] = mu;
	}
	_T[k-1][k-1] = (-1.0)*mu;
}





MixErDist::MixErDist(double m, double cv) : PhDist("MixEr", m, cv)
{
	if(cv >= 1){
		std::cerr << "ERROR (MixErDist::MixErDist): Invalid value for coefficient of variation" << std::endl;
		exit(1);
	}

	double cv2 = cv*cv;

	int k = (int)(1.0/cv2);

	double p = ( (k+1.0)/(1.0+cv2) ) 
					* ( cv2 - sqrt( (1.0-k*cv2)/(k+1.0) ) );

	double mu = ( p*k + (1.0-p)*(k+1) ) / _mean;

	// Initialize with a zero vector
	_alpha = RowVector(k+1);

	_alpha[0] = 1.0;
	
  // Initialize with a zero matrix
	_T = Matrix(k+1,k+1);

	for(int i=0; i < k-1; i++){
		_T[i][i] = (-1.0)*mu;
		_T[i][i+1] = mu;
	}
	_T[k-1][k-1] = (-1.0)*mu;
	_T[k-1][k] = (1.0-p)*mu;

	_T[k][k] = (-1.0)*mu;
}


BalHyperDist::BalHyperDist(double m, double cv) : PhDist("BalHyper", m, cv)
{
	if(cv <= 1){
		std::cerr << "ERROR (BalHyperDist::BalHyperDist): Invalid value for coefficient of variation" << std::endl;
		exit(1);
	}

	double cv2 = cv*cv;

  // Initialize with a zero vector
	_alpha = RowVector(2);
  // Initialize with a zero matrix
	_T = Matrix(2,2);

	double p = ( 1.0 + sqrt( (cv2-1.0)/(cv2+1.0) ) ) /2.0;

	double mu_1 = (2.0*p)/_mean;
	double mu_2 = (2.0*(1.0-p))/_mean;


	_alpha[0] = p;
	_alpha[1] = 1.0-p;

	_T[0][0] = (-1)*mu_1;
	_T[1][1] = (-1)*mu_2;
}


ColVector PhDist::unifVector(double tolerance, double gamma) const{

  if(gamma == 0.0){
    gamma = _T.getMaxAbsDiag();
  }

	double target_sum = gamma * _mean;

	Matrix I(_T.size(), _T.size()); 
	I.setIdentity();
	Matrix P = I + _T*(1.0/gamma);

	std::vector<double> col;

	ColVector e(_T.size()); 
	e.setAllOne();
	
	// alpha * P^n
	RowVector aPn = _alpha;

	double sum = 0;

	while(1){
		// alpha * P^n * e
		double prob = aPn * e;

		col.push_back(prob);
		sum += prob;

		if(target_sum - sum < tolerance){
			break;
		}else{
			aPn = aPn * P;
		}
	}

	return ColVector(col);
}

ColVector PhDist::unifDensityVector(double tolerance, double gamma) const{

  if(gamma == 0.0){
    gamma = _T.getMaxAbsDiag();
  }

	double target_sum = gamma * _mean;

	Matrix I(_T.size(), _T.size()); 
	I.setIdentity();
	Matrix P = I + _T*(1.0/gamma);

	std::vector<double> col;

	ColVector e(_T.size()); 
	e.setAllOne();

  ColVector T_e = (_T*(-1.0))*e;
	
	RowVector aPn = _alpha;

	double sum = 0;

	while(1){
		// alpha * P^n * e
		double prob = aPn * e;
    double density = aPn * T_e;

		col.push_back(density);
		sum += prob;

		if(target_sum - sum < tolerance){
			break;
		}else{
			aPn = aPn * P;
		}
	}

	return ColVector(col);
}



RowVector ConstDist::PoissonVec(double tolerance, double rate) const{
	return RowVector(PoissonFunc(rate*_mean, tolerance));
}

RowVector PhDist::PoissonVec(double tolerance, double rate) const{
	std::vector<double> row_vec;

	Matrix S = _T;
	RowVector beta = _alpha;

	Matrix I(S.size(), S.size());
	I.setIdentity();

	ColVector e(S.size());
	e.setAllOne();

	Matrix inv = (I - (1.0/rate)*S).getInverse();
	ColVector y = (1.0/rate) * inv * ((-1)*S)*e;


	row_vec.push_back(beta*y);
	double sum = row_vec[0];

	bool convergent = false;

	for(int i=1; row_vec[i-1] > 0; i++){
		y = inv*y;
		row_vec.push_back(beta*y);

		sum += row_vec[i];

		if(1.0 - sum < tolerance){
			convergent = true;
			break;
		}
	}

	if(!convergent){
		std::cerr << "ERROR (PoissonVec PH): Recursion not convergent" << std::endl;
		exit(1);
	}else{
		return RowVector(row_vec);
	}

}

RowVector ConstDist::PoissonVec_n(double tolerance, int n, double rate) const{
	return RowVector(PoissonFunc_n(rate*_mean, tolerance, n));
}

RowVector PhDist::PoissonVec_n(int n, double rate) const
{
	std::vector<double> row_vec;

	Matrix S = _T;
	RowVector beta = _alpha;

	Matrix I(S.size(), S.size());
	I.setIdentity();

	ColVector e(S.size());
	e.setAllOne();

	Matrix inv = (I - (1.0/rate)*S).getInverse();
	ColVector y = (1.0/rate) * inv * ((-1)*S)*e;


	row_vec.push_back(beta*y);
	for(int i=1; i < n; i++){
		y = inv*y;
		row_vec.push_back(beta*y);

	}

	return RowVector(row_vec);
}

std::vector<double> PoissonFunc(double lmdt, double tolerance){
	int base = lmdt;

	std::vector<double> undist_lower;
	std::vector<double> undist_upper;

	double sum_lower = 0;
	double sum_upper = 0;

	double value_lower = 1;
	double value_upper = 1;

	bool stop_lower = false;
	bool stop_upper = false;

	bool isConvergent = false;

	for(int i = 0; value_lower > 0 || value_upper > 0; i++){
		if(base - i <= 0){stop_lower = true;}
		if(!stop_lower){
			double r = (base-i)/lmdt;

			// Utilize the fact that value_lower * (1-r^k)/(1-r) < torelance => value_lower * 1/(1-r) < torelance 
			if(value_lower < tolerance * (1.0-r)){
				stop_lower = true;
			}

			value_lower = value_lower * r; 
			sum_lower = sum_lower + value_lower;

			undist_lower.push_back(value_lower);
		}
		if(!stop_upper){
			double r = lmdt/(base+i+1);

			// Utilize the fact that value_lower * 1/(1-r) < torelance 
			if(value_upper < tolerance * (1.0-r)){
				stop_upper = true;

			}

			value_upper = value_upper * r;
			sum_upper = sum_upper + value_upper;

			undist_upper.push_back(value_upper);
		}

		if(stop_lower && stop_upper){
			isConvergent = true;
			break;
		}
	}
	if(isConvergent){
		while((int)undist_lower.size() < base){
			undist_lower.push_back(0.0);
		}

		double sum = 1 + sum_lower + sum_upper;
		std::vector<double> dist;

		for(int k = undist_lower.size()-1; k >= 0;  k--){
			dist.push_back(undist_lower[k]/sum);
		}

		dist.push_back(1.0/sum);

		for(unsigned int k = 0; k < undist_upper.size(); k++){
			dist.push_back(undist_upper[k]/sum);
		}

		return dist;
	} else{
		std::cerr << "ERROR (PoissonFunc): Recursion is not convergent" << std::endl;
		exit(1);
	}
}

std::vector<double> PoissonFunc_n(double lmdt, double tolerance, int n)
{
	std::vector<double> dist 
		= PoissonFunc(lmdt, tolerance);

	std::vector<double> dist_n;
	for(int i=0; i < n && (unsigned int)i < dist.size(); i++){
		dist_n.push_back(dist[i]); 
	}

	return dist_n;
}

std::vector<double> BinomialFunc(int k, double p, double tolerance, int& start, int& end){
	std::vector<long double> dist(k+1, 0.0);

	int ave = k*p;

	long double log_p_ave = 0;

	if(p > 0.5){
		for(int i = 1; i <= k-ave; i++){
			log_p_ave = log_p_ave + std::logl( (long double)(1.0-p) * (long double)(ave+i) / i);
		}
		log_p_ave = log_p_ave + ave*std::logl(p);
	} else{
		for(int i = 0; i < ave; i++){
			log_p_ave = log_p_ave + std::logl( (long double)p * (long double)(k-i) / (ave-i));
		}
		log_p_ave = log_p_ave + (k-ave)*std::logl(1.0-p);
	}

	
	dist[ave] = std::expl(log_p_ave);

	int i_fwd = ave;
	long double value_fwd = dist[ave];

	int i_bak = ave;
	long double value_bak = dist[ave];

	bool fwd_stop = (i_fwd == k);
	bool bak_stop = (i_bak == 0);

	bool convergent = false;

	long double sum = dist[ave];

	while(1){
		bool calcfwd = !fwd_stop && (bak_stop || (value_fwd > value_bak) );

		if(calcfwd){
			value_fwd = dist[i_fwd] 
							* (p/(1.0-p)) * ((long double)(k-i_fwd)/(i_fwd+1));
			dist[i_fwd+1] = value_fwd;
			sum = sum + value_fwd;
			i_fwd++;
			if(i_fwd == k){fwd_stop = true;}
		}else{
			value_bak = dist[i_bak]
							* ((1.0-p)/p) * ((long double)i_bak/(k-i_bak+1));
			dist[i_bak-1] = value_bak;
			sum = sum + value_bak;

			i_bak--;				
			if(i_bak == 0){bak_stop = true;}
		}

		if(1.0 - sum < tolerance){convergent = true; break;}

		if(bak_stop && fwd_stop){break;}
	}

	if(!convergent){
		std::cerr << "ERROR (BinomialFunc(" << k << "," << p
			<< ")): Recursion is not convergent" << std::endl;
		exit(1);
	}
	
	start = i_bak;
	end = i_fwd;

  std::vector<double> dist_double;
  for(int i=0; i < dist.size(); i++){
    dist_double.push_back( (double)dist[i] );
  }

	return dist_double;
}

std::vector<double> BinomialFuncFull(int k, double p){
	std::vector<double> dist(k+1, 0.0);

	int ave = k*p;

	double log_p_ave = 0;

	if(p > 0.5){
		for(int i = 1; i <= k-ave; i++){
			log_p_ave = log_p_ave + std::logl( (1.0-p) * (double)(ave+i) / i);
		}
		log_p_ave = log_p_ave + ave*std::logl(p);
	} else{
		for(int i = 0; i < ave; i++){
			log_p_ave = log_p_ave + std::logl( p * (double)(k-i) / (ave-i));
		}
		log_p_ave = log_p_ave + (k-ave)*std::logl(1.0-p);
	}

	
	dist[ave] = std::expl(log_p_ave);

	int i_fwd = ave;
	double value_fwd = dist[ave];

	int i_bak = ave;
	double value_bak = dist[ave];

	bool fwd_stop = (i_fwd == k);
	bool bak_stop = (i_bak == 0);

	while(true){
    fwd_stop = (i_fwd == k);
	  bak_stop = (i_bak == 0);

		if(bak_stop && fwd_stop){break;}

		bool calcfwd = !fwd_stop && (bak_stop || (value_fwd > value_bak) );

		if(calcfwd){
			value_fwd = dist[i_fwd] 
							* (p/(1.0-p)) * ((double)(k-i_fwd)/(i_fwd+1));
			dist[i_fwd+1] = value_fwd;
			i_fwd++;
		}else{
			value_bak = dist[i_bak]
							* ((1.0-p)/p) * ((double)i_bak/(k-i_bak+1));
			dist[i_bak-1] = value_bak;
			i_bak--;				
		}

	}
	return dist;
}




RowVector Distribution::PoissonVec(double tolerance, double rate) const
{
	std::cerr << "ERROR (Distribution::PoissonVec(double tolerance, double rate) : Computational method is undefined for this type of distribution" << std::endl;
	exit(1);
}

RowVector Distribution::PoissonVec(double tolerance, double rate, bool show_err) const
{
	std::cerr << "ERROR (Distribution::PoissonVec(double tolerance, double rate, bool show_err)) : Computational method is undefined for this type of distribution" << std::endl;
	exit(1);
}

RowVector Distribution::PoissonVec_n(double tolerance, int n, double rate) const 
{
	std::cerr << "ERROR (Distribution::PoissonVec_n(double tolerance, int n, double rate)) : Computational method is undefined for this type of distribution" << std::endl;
	exit(1);
}

RowVector Distribution::PoissonVec_n(int n, double rate) const
{
	std::cerr << "ERROR (PoissonVec_n(int n, double rate)) : Computational method is undefined for this type of distribution" << std::endl;
	exit(1);
}

RowVector Distribution::PoissonVec_n(int n, double rate, bool show_err, double& err) const
{
	std::cerr << "ERROR (PoissonVec_n(int n, double rate, bool show_err)) : Computational method is undefined for this type of distribution" << std::endl;
	exit(1);
}
