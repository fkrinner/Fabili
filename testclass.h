#ifndef TESTCLASS
#define TESTCLASS
#include"fabiliFunction.h"

class testclass : public fabiliFunction {
	public:
		testclass(double c, std::vector<double> b, std::vector<std::vector<double> > a);

		double scalarEval(const std::vector<double>& parameters) const override;
		evalType eval(const std::vector<double>& parameters) const override;

		size_t dim() const override {return _dim;}
	protected:
		size_t                            _dim;
		double                            _C;
		std::vector<double>               _B;
		std::vector<std::vector<double> > _A;
};

#endif//TESTCLASS
