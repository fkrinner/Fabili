#ifndef TESTDERIVATIVES
#define TESTDERIVATIVES
#include<memory>
#include<vector>
#include"fabiliFunction.h"

void testDerivatives(std::shared_ptr<fabiliFunction> function, std::vector<double> point,  double delta = 1.e-5, bool testHessian = true) {
	size_t dim = function->dim();
	if (dim != point.size()) {
		std::cout << "Number of parameters does not match" << std::endl;
		throw;
	}
	evalType eval = function->eval(point);
	for (size_t i = 0; i < dim; ++i) {
		point[i] += delta;
		evalType delEval = function->eval(point);
		point[i] -= delta;
		std::cout << "gradient(" << i << "): coded: " << eval.gradient(i) << "; numerical: " << (delEval.value - eval.value)/delta << std::endl;
		if (testHessian) {
			for (size_t j = 0; j < dim; ++j) {
				std::cout << "hessian(" << i << "," << j << "): coded: " << eval.hessian(i,j) << "; numerical: " << (delEval.gradient(j) - eval.gradient(j))/delta << std::endl;
			}
		}
	}
}
#endif//TESTDERIVATIVES
