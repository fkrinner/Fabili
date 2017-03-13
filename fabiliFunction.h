#ifndef FABILI_FUNCTION
#define FABILI_FUNCTION
#include<vector>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>

struct evalType {
	double          value;
	Eigen::VectorXd gradient;
	Eigen::MatrixXd hessian;
};


class fabiliFunction {
	public:
		fabiliFunction(): _dim(0){};
		virtual evalType eval(const std::vector<double>& parameters) const = 0;

		size_t dim() const {return _dim;}
	protected:
		size_t _dim;
};
#endif//FABILI_FUNCTION
