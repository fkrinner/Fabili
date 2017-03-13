#include"fabili.h"

fabili::fabili(std::shared_ptr<fabiliFunction> function) : _numericalLimit(1.e-8), _stepSize(0.), _function(function) {};

std::vector<double> fabili::estimate(const std::vector<double>& point) const {
	const size_t dim = point.size();

	Eigen::VectorXd P(dim);
	for (size_t i = 0; i < dim; ++i) {
		P(i) = point[i];
	}

	evalType evalPoint = _function->eval(point);

	Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(evalPoint.hessian);
	Eigen::MatrixXcd V = eigenSolver.eigenvectors().transpose();

	Eigen::VectorXcd Drot = V*evalPoint.gradient;
	Eigen::VectorXcd Prot = V*P;
	Eigen::VectorXd estimateRot(dim);

	for (size_t i = 0; i < dim; ++i) {
		estimateRot(i) = estimate1D(eigenSolver.eigenvalues()[i], Drot(i), Prot(i));
	}
	Eigen::VectorXcd estimate = V.adjoint() * estimateRot;

	std::vector<double> retVal(dim);
	for(size_t i = 0; i < dim; ++i) {
		std::complex<double> val = estimate(i);
		if (not isZero(val.imag())) {
			std::cout << "fabili::estimate(...): ERROR: Complex estimate encountered... using old value" << std::endl;
			retVal[i] = point[i];
		} else {
			retVal[i] = val.real();
		}
	}
	return retVal;	
}

double fabili::estimate1D(const std::complex<double>& DD, const std::complex<double>& D, const std::complex<double>& x) const {
		if (not isZero(DD.imag())) { 
			std::cout << "fabili::esitmate1D(...): ERROR: Second derivative is not real " << DD << ". Return old value" << std::endl;
			return x.real();
		}
		double DDf = DD.real();
		if (not isZero(D.imag())) {
			std::cout << "fabili::esitmate1D(...): ERROR: First derivative is not real " << D << ". Return old value" << std::endl;
			return x.real();
		}
		double Df = D.real();
		if (not isZero(x.imag())) {
			std::cout << "fabili::esitmate1D(...): ERROR: Point is not real " << x << ". Return real part" << std::endl;
			return x.real();
		}
		double X   = x.real();
		double A   = DDf/2;
		double B   = Df - DDf * X;
		if (isZero(DDf)) {
			return handleZeroEigenvalue(A, B, X);
		} else if (DDf > 0.) {
			return handlePositiveEigenvalue(A, B, X);
		}
		return handleNegativeEigenvalue(A, B, X);
}

double fabili::handleZeroEigenvalue(double A, double B, double X) const {
	if (isZero(B)) {
		return X;
	} else if (B < 0.) {
		return X + _stepSize;
	}
	return X - _stepSize;
}

double fabili::handlePositiveEigenvalue(double A, double B, double X) const {
	return -B/A/2.;
}

double fabili::handleNegativeEigenvalue(double A, double B, double X) const {
	double BB = 4*A*X + B;
	return BB/(2*A); // AA = -A
}

bool fabili::isZero(double val) const {
	if (abs(val) < _numericalLimit) {
		return true;
	}
	return false;
}



