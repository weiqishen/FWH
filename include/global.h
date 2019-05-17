// This is a header file which includes all the global variables and all the functions.
#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "ndarray.h"

#define VERSION 1.0

struct Params
{
	double T_static,p_static;
	double gamma;
	double R_gas;
	string fwh_surf_fname;
	string output_fname;
	int endcap_avg;
	size_t n_endcaps;
	ndarray<double> endcap_x;
};

struct FWH_surf
{
	//data
	ndarray<double> data;

	//gemetry
	ndarray<double> A;
	ndarray<double> normal;
	ndarray<double> center;
	size_t n_eles;

	double dt;
	ndarray<double> tau; //source time
};

struct Observer
{
	//observer position
	ndarray<double> x; 
	ndarray<double> y; 
	ndarray<double> z; 
	size_t n_oberver;
};

struct dist
{
	ndarray<double> r;		 //distance vector
	ndarray<double> mag_r; //distance magnitude
};

extern double PI;
extern Params input;
extern FWH_surf faces;//declare fwh surface
extern Observer microphone;//declare observation microphone

