#include"fabili.h"
#include"testclass.h"
#include"fabiliFunction.h"
#include"testDerivatives.h"
#include<memory>
#include<vector>

int main() {
	double C                            = -123.;
	std::vector<double> B               = {12., -321., .01};
	std::vector<std::vector<double> > A = {	{ 12., 0.1, 0.12},
						{ -.2, 18., -.1 },
						{ -.1, 0.12, 100.}};
	
	B = {-123.,0.,0.};

	A = {{30.,-1.,0.,},{-1.,20.,0.},{0.,0.,10.}};

	const size_t dim = B.size();
	std::shared_ptr<testclass> func = std::make_shared<testclass>(C, B, A);
	fabili fab(func);

	std::vector<double> point = {2.,2.,1.};
	std::vector<double> estim = fab.minimize(point, 100).second;
	evalType atPoint = func->eval(point);
	evalType atEstim = func->eval(estim);

//	testDerivatives(func, point);
	return 0;
}
