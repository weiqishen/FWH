#pragma once

#include "global.h"

class F1A_solver
{
public:
	F1A_solver();
	~F1A_solver();
	void solve(void);

private:
	//funtion member
	void calc_dis_src_obs(size_t obs_id);
	double get_begin_signal_time();
	double get_final_signal_time();
	void calc_pressure_term(double endt, double begint);
	double get_weight(double x);
	void write(size_t obs_id);
	//data member
	double c_infty,rho_infty;								 //speed of sound
	ndarray<double> t;		 //observer time
	ndarray<double> pTotal;  // Total pressure at the microphone
	dist src2ob;				//vector of sources to observer
};
