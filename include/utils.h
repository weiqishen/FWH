inline double linearInterpolate(double tL, double tU, double t, double pL, double pU)
{
		return  pL + (t - tL) * (pU - pL) / (tU - tL);
}
inline double calc_slope(double pL, double pU, double h)
{
	return (pU - pL) / h;
}
