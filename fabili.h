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

		std::vector<double> estimate(const std::vector<double>& point) const;

		double estimate1D(const std::complex<double>& DD, const std::complex<double>& D, const std::complex<double>& x) const;
		double handleZeroEigenvalue(double A, double B, double X) const;
		double handlePositiveEigenvalue(double A, double B, double X) const;
		double handleNegativeEigenvalue(double A, double B, double X) const;

		bool isZero(double val) const;

	protected:
		double                          _numericalLimit;
		double                          _stepSize;
		std::shared_ptr<fabiliFunction> _function;

};

#endif//FABILI__ILIBAF
