#include"fabili.h"

fabili::fabili(std::shared_ptr<fabiliFunction> function) : _dim(function->dim()), _numericalLimit(1.e-8), _stepSize(0.), _function(function) {};

std::pair<bool, std::vector<double> > fabili::minimize(const std::vector<double>& startingPoint, int callLimit) const {
	int count = 0;
	std::vector<double> point(startingPoint);
	while (true) {
		std::pair<bool, std::vector<double> > estim = estimate(point);
		if (estim.first) {
			std::cout << "Minimization finished after " << count << " function calls." << std::endl;
			return estim;
		}
		++count;
		point = estim.second;
		if (count > callLimit && callLimit > 0) {
			break;
		}
	}
	std::cout << "Minimization unsuccessful, call limit of " << callLimit << " function calls exceeded" << std::endl;
	return std::pair<bool, std::vector<double> >(false, point);
} 

std::pair<bool, std::vector<double> > fabili::estimate(const std::vector<double>& point) const {
	const size_t dim = point.size();

	Eigen::VectorXd P(dim);
	for (size_t i = 0; i < dim; ++i) {
		P(i) = point[i];
	}

	evalType evalPoint = _function->eval(point);

	if (stoppingCriterion(evalPoint)) {
		return std::pair<bool, std::vector<double> >(true, point);
	}

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
	return std::pair<bool, std::vector<double> >(false, retVal);	
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
			std::cout << "Zero handler for " << A << std::endl;
			return handleZeroEigenvalue(A, B, X);
		} else if (DDf > 0.) {
			return handlePositiveEigenvalue(A, B, X);
		}
		return handleNegativeEigenvalue(A, B, X);
}

double fabili::handleZeroEigenvalue(double A, double B, double X) const {
std::cout << "Zero handler called" << std::endl;
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
std::cout << "Negative handler called" << std::endl;
	double BB = 4*A*X + B;
// f(x) = Ax^2 + Bx + C with A < 0
// Predict for negative curvature, as for a parabula with same gradient but opposite sign curvature:
// g(x) = f(x), g'(x) = f'(x), g''(x) = -f''(x) for x = X
// g(x) = A'x^2 + B'x + C' => Prediction = -B'/(2A') 
// A' = -A; B' = 4AX + B; C' is irrelevant for the prediction
	return BB/(2*A);
}

bool fabili::isZero(double val) const {
	if (fabs(val) < _numericalLimit) {
		return true;
	}
	return false;
}

bool fabili::stoppingCriterion(const evalType& point) const {
	double absGrad = 0.;
	for (size_t i = 0; i < _dim; ++i) {
		absGrad += pow(point.gradient(i),2);
	}
	absGrad = pow(absGrad,.5);
	std::cout << "|f'(x)| = " << absGrad << " (f(x) = " << point.value << ")" << std::endl;
	if (absGrad < _numericalLimit) {
		return true;
	}
	return false;
}

