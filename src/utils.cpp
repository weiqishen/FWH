#include "utils.h"

using namespace std;

double linearInterpolate(double tL, double tU, double t, double pL, double pU)
{

	double p;

	if (t == 0)
		p = 0;
	else
		p = pL + ((t - tL) * (pU - pL) / (tU - tL));

	return p;
}

double secondOrderCentralDiff(double pL, double pU, double h)
{
	double result;

	result = (pU - pL) / (2 * h);

	return result;
}

double secondOrderBackwardDiff(double p, double pL, double pLL, double h)
{
	double result;

	result = (((3.0 / 2.0) * p) - (2.0 * pL) + (0.5 * pLL)) / h;

	return result;
}