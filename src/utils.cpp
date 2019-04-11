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

double get_weight(double x)
{
	if (x < 25 * 0.0508)
	{
		return 1.0;
	}
	else if (x == 25 * 0.0508 || x == 26 * 0.0508 || x == 27 * 0.0508 || x == 28 * 0.0508 || x == 29 * 0.0508 || x == 30 * 0.0508)
	{
		return 1. / 6.;
	}
	else if (x > 25 * 0.0508 && x < 26 * 0.0508)
	{
		return 5. / 6.;
	}
	else if (x > 26 * 0.0508 && x < 27 * 0.0508)
	{
		return 4. / 6.;
	}
	else if (x > 27 * 0.0508 && x < 28 * 0.0508)
	{
		return 3. / 6.;
	}
	else if (x > 28 * 0.0508 && x < 29 * 0.0508)
	{
		return 2. / 6.;
	}
	else if (x > 29 * 0.0508 && x < 30 * 0.0508)
	{
		return 1. / 6.;
	}
	else
	{
		return 0.;
	}
}