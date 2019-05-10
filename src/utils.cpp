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

double calc_slope(double pL, double pU, double h)
{
	double result;

	result = (pU - pL) / (2 * h);

	return result;
}