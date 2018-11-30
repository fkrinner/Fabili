#ifndef FABILI__ILIBAF
#define FABILI__ILIBAF

#include<iostream>
#include<vector>
#include<memory>

#include"fabiliFunction.h"

#include<Eigen/Dense>
#include<Eigen/Eigenvalues>

class fabili {
	public:
		fabili(std::shared_ptr<fabiliFunction> function);

		std::pair<bool, std::vector<double> > estimate(const std::vector<double>& point) const;

		std::pair<bool, std::vector<double> > minimize(const std::vector<double>& startingPoint, size_t callLimit = 0) const; 

		double estimate1D(const std::complex<double>& DD, const std::complex<double>& D, const std::complex<double>& x) const;
		double handleZeroEigenvalue(double A, double B, double X)     const;
		double handlePositiveEigenvalue(double A, double B, double X) const;
		double handleNegativeEigenvalue(double A, double B, double X) const;
		std::vector<double> handleWorsePoint(const std::vector<double>& point, const evalType& evalPoint) const;

		bool isZero(double val) const;
		bool stoppingCriterion(const evalType& point) const;
		bool setStepSize(double value) {_stepSize = value; return true;}
		bool setNstepWorse(size_t value) {_nStepWorse = value; return true;}

	protected:
		size_t                          _dim;
		size_t                          _nStepWorse; // If estimated point is worse, scan the area in between with _nStepWorse steps
		size_t                          _maxCall;
		mutable size_t                  _nCall;
		double                          _numericalLimit;
		double                          _stepSize;
		mutable double                  _lastEval;
		mutable std::vector<double>     _lastParameters;
		std::shared_ptr<fabiliFunction> _function;

};

#endif//FABILI__ILIBAF
