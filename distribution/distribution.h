#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include<cmath>
#include<iostream>
#include<vector>
#include<string>

#include"../matrix/matrix.h"
#include"../matrix/myVector.h"


class Distribution{
	public:
		// Type of the distribution
		std::string _type;

		// Mean
		double _mean;
		// Coefficient of variation
		double _coef_v;

		Distribution(){};

		Distribution(std::string tp, double m, double cv);
		Distribution(const char* tp, double m, double cv);

		virtual ~Distribution(){};

		// Returns the probability mass function of the number of Poisson arrivals in the random interval
		virtual RowVector PoissonVec(double tolerance, double rate, bool show_err) const;
		virtual RowVector PoissonVec(double tolerance, double rate) const;
		virtual RowVector PoissonVec_n(double tolerance, int n, double rate) const;
		virtual RowVector PoissonVec_n(int n, double rate, bool show_err, double& err) const;
		virtual RowVector PoissonVec_n(int n, double rate) const;

		virtual void print() const;
};

// Deterministic distribution
class ConstDist : public Distribution{
	public:
		ConstDist(double m) : Distribution("Const", m, 0.0){}

		RowVector PoissonVec(double tolerance, double rate) const;
		RowVector PoissonVec_n(double tolerance, int n, double rate) const;
};

// Phase-type distribution
class PhDist : public Distribution{
	public:
		RowVector _alpha;
		Matrix _T;

		PhDist(RowVector vec, Matrix mat);
		PhDist(std::string tp, double m, double cv)
									: Distribution( tp, m, cv){}
		PhDist(const char* tp, double m, double cv)
									: Distribution( tp, m, cv){}

		// Returns a uniformization vector
		ColVector unifVector(double tolerance, double gamma=0.0) const;
    ColVector unifDensityVector(double tolerance, double gamma=0.0) const;

		RowVector PoissonVec(double tolerance, double rate) const;
		RowVector PoissonVec_n(int n, double rate) const;
};

// Exponential distribution
class ExpDist : public PhDist{
	public:
    // m: Mean
		ExpDist(double m) : PhDist("Exp", m, 1.0){
			_alpha = RowVector(1);
			_alpha[0] = 1.0;

			_T = Matrix(1,1);
			_T[0][0] = -1.0/_mean;
		}
};

// Erlang distribution
class ErDist : public PhDist{
	public:
		int _stage;
		// m: Mean, k: number of stages
		ErDist(double m, int k);

		void print() const;
};


// Mixture of k-stage and (k+1)-stage Erlang distributions
class MixErDist : public PhDist{
	public:
		// m: mean, cv: coefficient of variation
		MixErDist(double m, double cv);
};


// Balanced hyper-exponential distribution
class BalHyperDist : public PhDist{
	public:
		// m: mean, cv: coefficient of variation
		BalHyperDist(double m, double cv);
};


// Probability mass function of Poisson distribution with mean lmdt
// (Truncation point is determined based on the tolerance)
std::vector<double> PoissonFunc(double lmdt, double tolerance);
// (Truncation point is specified by n)
std::vector<double> PoissonFunc_n(double lmdt, double tolerance, int n);

// Binomial distribution with k trials and success probability p (with trunction)
// The range of non-zero elements are substituted into start and end
std::vector<double> BinomialFunc(int k, double p, double tolerance,
													int& start, int& end);

// Binomial distribution with k trials and success probability p (without truncation)
std::vector<double> BinomialFuncFull(int k, double p);

#endif
