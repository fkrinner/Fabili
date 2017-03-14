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
		fabiliFunction(){};
		virtual evalType eval(const std::vector<double>& parameters) const = 0;

		virtual size_t dim() const {return 0;}
};
#endif//FABILI_FUNCTION
