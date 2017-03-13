#ifndef TESTCLASS
#define TESTCLASS
#include"fabiliFunction.h"

class testclass : public fabiliFunction {
	public:
		testclass(double c, std::vector<double> b, std::vector<std::vector<double> > a);

		evalType eval(const std::vector<double>& parameters) const;
	protected:
		double                            _C;
		std::vector<double>               _B;
		std::vector<std::vector<double> > _A;
};

#endif//TESTCLASS
