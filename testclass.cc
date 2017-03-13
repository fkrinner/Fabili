#include"testclass.h"
#include<iostream>
testclass::testclass(double c, std::vector<double> b, std::vector<std::vector<double> > a) : fabiliFunction(), _C(c), _B(b), _A(a) {
	_dim = _B.size();
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
	evalType retVal;
	retVal.value     = _C;
	retVal.gradient = Eigen::VectorXd::Zero(_dim);
	retVal.hessian   = Eigen::MatrixXd::Zero(_dim, _dim);
	for (size_t i = 0; i < _dim; ++i) {
		retVal.value += _B[i] * parameters[i];
		retVal.gradient(i) += _B[i];
		for (size_t j = 0; j < _dim; ++j) {
			retVal.value += parameters[i] *  _A[i][j] * parameters[j];
			retVal.gradient(i)  += (_A[i][j] + _A[j][i]) * parameters[j];
			retVal.hessian(i,j) += (_A[i][j] + _A[j][i]);
		}
	}
	return retVal;
}
