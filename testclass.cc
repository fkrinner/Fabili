#include"testclass.h"
#include<iostream>

double f(double x) {
	return pow(x,2);
}

double Df(double x) {
	return 2*x;
}

double DDf(double x) {
	return 2.;
}

testclass::testclass(double c, std::vector<double> b, std::vector<std::vector<double> > a) : fabiliFunction(), _dim(b.size()), _C(c), _B(b), _A(a) {
	if (_A.size() != _B.size()) {
		std::cout << "A and B sizes do not match" << std::endl;
		throw;
	}
	for (size_t i = 0; i < _dim; ++i) {
		if (_A[i].size() != _dim) {
			std::cout << "line " << i << " of A has the wrong size" << std::endl;
			throw;
		}
	}
}

evalType testclass::eval(const std::vector<double>& parameters) const {
	if (parameters.size() != _dim) {
		std::cout << "Wrong size of parameters" << std::endl;
		throw;
	}
	double X = _C; // X = xAx + Bx + C
	std::vector<double> baseGradient(_dim,0.); // dX/dx_i
	for (size_t i = 0; i < _dim; ++i) {
		X += _B[i] * parameters[i];
		baseGradient[i] += _B[i];
		for (size_t j = 0; j < _dim; ++j) {
			baseGradient[i] += (_A[i][j] + _A[j][i]) * parameters[j];
			X               += parameters[i] *  _A[i][j] * parameters[j];
		}
	} 

	evalType retVal;
	retVal.value        = f(X); // f(X)
	retVal.gradient     = Eigen::VectorXd::Zero(_dim);
	retVal.hessian      = Eigen::MatrixXd::Zero(_dim, _dim);

	for (size_t i = 0; i < _dim; ++i) {
		retVal.gradient(i) += Df(X)*baseGradient[i]; // f'(X) dX/dx_i
		for (size_t j = 0; j < _dim; ++j) {
			retVal.hessian(i,j) += Df(X)*(_A[i][j] + _A[j][i]); // f'(X) d2X/dx_idx_j
			retVal.hessian(i,j) += DDf(X)*baseGradient[i]*baseGradient[j]; // f''(X) dX/dx_i dX/dx_j
		}
	}
	return retVal;
}
