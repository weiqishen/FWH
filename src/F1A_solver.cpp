// This C program calculates the SPL per Hz. at the observer position from the source position.

#include "global.h"
#include "F1A_solver.h"
#include "utils.h"
#include <numeric>
#include <fstream>
#include <sys/stat.h>

using namespace std;
F1A_solver::F1A_solver()
{
	c = sqrt(input.gamma * input.R_gas * input.T_static);
}

F1A_solver::~F1A_solver()
{
}

void F1A_solver::solve(void)
{
	double begint;
	double endt;
	struct stat ex;

	//create folder if not exist
	if (stat(input.output_fname.c_str(), &ex) == -1)
		mkdir(input.output_fname.c_str(), 0755);

	//for each observer
	for (size_t i = 0; i < microphone.n_oberver; i++)
	{
		cout << "Calculating the terms of pressure at observer location " << i << " of " << microphone.n_oberver << " as given in Farassat's 1A formulation..." << endl;

		// Calculating the distance between the source and the observer points.
		calc_dis_src_obs(i);
		// Calculating when the first sound signal will reach the farthest microphone/s.
		begint = get_begin_signal_time();
		// Calculating when the final sound signal will reach the farthest microphone/s.
		endt = get_final_signal_time();

		if (endt - begint < 0)
		{
			cout << "  -> ****Initial pressure signal from the farthest face reaches the observer after all the pressure signals from the nearest face have reached there." << endl;
			cout << "  -> Begin time: " << begint << "\t"
				 << "End time: " << endt << endl;
			cout << "  -> ****Skipping the calculation for this microphone location" << i << "." << endl;
			continue;
		}

		// Calculating the source time based on the observer time and interpolating the pressure values on the surface.
		//calculate pressure term in F1A formulation
		calc_pressure_term(endt, begint);

		//write file
		write(i);
	}
	cout << "Done." << endl;
}

void F1A_solver::calc_dis_src_obs(size_t obs_id)
{
	cout << " -> Calculating the distance between source and microphone " << obs_id << " ... " << flush;

	// Allocates the memory to store the distance between the microphone and the source faces (part 2).
	src2ob.r.setup({3, faces.n_eles});
	src2ob.mag_r.setup(faces.n_eles);

	// Calculates the difference between the co-ordinates between microphone and source faces
	// and also finds the magnitude of differnece using that.
	for (size_t j = 0; j < faces.n_eles; j++)
	{
		src2ob.r({0, j}) = microphone.x(obs_id) - faces.center({0, j});
		src2ob.r({1, j}) = microphone.y(obs_id) - faces.center({1, j});
		src2ob.r({2, j}) = microphone.z(obs_id) - faces.center({2, j});
		src2ob.mag_r(j) = sqrt(inner_product(src2ob.r.get_ptr({0, j}), src2ob.r.get_ptr({0, j + 1}), src2ob.r.get_ptr({0, j}), 0.));
	}

	cout << " Done.\n";
}

double F1A_solver::get_begin_signal_time()
{
	return faces.tau(0) + src2ob.mag_r.get_max() / c;
}

double F1A_solver::get_final_signal_time()
{
	return faces.tau(faces.tau.get_len() - 1) + src2ob.mag_r.get_min() / c;
}

void F1A_solver::calc_pressure_term(double endt, double begint)
{
	//element local arrays
	ndarray<double> Lr;		 // Lr = Li*\hat{r}/|r| -> Li = p\hat{n} and r = radiation vector from source faces to microphone location
	ndarray<double> rhoCalc; // Interpolated density at source
	ndarray<double> Machr;   // M*\hat{r}/|r| of flow -> M = Mach number and r = radiation vector from source faces to mic location
	ndarray<double> un;		 // u*\hat{n} of flow -> Normal velocity
	ndarray<double> ur;		 // u*\hat{r}/|r| of flow -> velocity in radiation direction

	//temporary variables
	size_t counter;//global time counters
	double tauCalc, pressure;//interpolate source time and source pressure
	ndarray<double> L(3), velocity(3);
	double dlr_dt;
	double dMachr_dt;
	double dun_dt;
	double drho_dt;
	double w, p1, p2, p3, p4;
	bool flag;

	cout << " -> Calculating variables at interpolated source time... " << endl;

	//initialize global arrays
	//observer time array
	t.setup((size_t)((endt - begint) / (0.5 * faces.dt)) + 1);
	for (size_t i = 0; i < t.get_len(); i++)
		t(i) = begint + i * 0.5 * faces.dt;

	//final pressure array
	pTotal.setup(t.get_len());
	pTotal = 0.;

	//initialize element local arrays
	Lr.setup(t.get_len());
	rhoCalc.setup(t.get_len());
	Machr.setup(t.get_len());
	un.setup(t.get_len());
	ur.setup(t.get_len());

	//for each source
	for (size_t j = 0; j < faces.n_eles; j++)
	{
		cout << setprecision(3) << (double)j / faces.n_eles * 100 << "\%   \r" << flush;
		
		counter = 0;
		//calculate interpolated source time variables
		for (size_t timeloop = 0; timeloop < t.get_len(); timeloop++) //for each observer time step
		{
			tauCalc = t(timeloop) - src2ob.mag_r(j) / c; //the interpolated src time

			//calculate interpolated variables
			flag = false;
			while (counter < faces.tau.get_len() - 1)
			{
				if (tauCalc >= faces.tau(counter) && tauCalc <= faces.tau(counter + 1))
				{
					pressure = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, faces.data({counter, j, 4}), faces.data({counter + 1, j, 4}));
					//p\hat{n}
					L(0) = pressure * faces.normal({0, j});
					L(1) = pressure * faces.normal({1, j});
					L(2) = pressure * faces.normal({2, j});
					//p\hat{n}*\hat{r}/|r|, project normal stress on distance dir
					Lr(timeloop) = inner_product(L.get_ptr(), L.get_ptr(3), src2ob.r.get_ptr({0, j}), 0.) / src2ob.mag_r(j);

					rhoCalc(timeloop) = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, faces.data({counter, j, 0}), faces.data({counter + 1, j, 0}));

					velocity(0) = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, faces.data({counter, j, 1}), faces.data({counter + 1, j, 1}));
					velocity(1) = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, faces.data({counter, j, 2}), faces.data({counter + 1, j, 2}));
					velocity(2) = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, faces.data({counter, j, 3}), faces.data({counter + 1, j, 3}));

					//\hat{u}/c*\hat{r}/|r| project mach number on distance dir
					Machr(timeloop) = inner_product(velocity.get_ptr(), velocity.get_ptr(3), src2ob.r.get_ptr({0, j}), 0.) / src2ob.mag_r(j) / c;
					//\hat{u}*\hat{n} normal velocity
					un(timeloop) = inner_product(velocity.get_ptr(), velocity.get_ptr(3), faces.normal.get_ptr({0, j}), 0.);
					//|Mr|*c project velocity on distance dir
					ur(timeloop) = Machr(timeloop) * c;
					flag=true;
					break;
				}
				counter++;
			}
			if (flag == false)
				Fatal_Error("cant find interpolated source time");
		}

		//add source term of each element to each observer time step
		for (size_t timeloop = 0; timeloop < t.get_len(); timeloop++) //for each observer time step
		{
			if (timeloop > 0 && timeloop < t.get_len() - 1)
			{
				dlr_dt = secondOrderCentralDiff(Lr(timeloop - 1), Lr(timeloop + 1), 0.5 * faces.dt);
				dMachr_dt = secondOrderCentralDiff(Machr(timeloop - 1), Machr(timeloop + 1), 0.5 * faces.dt);
				dun_dt = secondOrderCentralDiff(un(timeloop - 1), un(timeloop + 1), 0.5 * faces.dt);
				drho_dt = secondOrderCentralDiff(rhoCalc(timeloop - 1), rhoCalc(timeloop + 1), 0.5 * faces.dt);
			}
			else if (timeloop == 0)
			{
				dlr_dt = secondOrderBackwardDiff(-Lr(timeloop), -Lr(timeloop + 1), -Lr(timeloop + 2), 0.5 * faces.dt);
				dMachr_dt = secondOrderBackwardDiff(-Machr(timeloop), -Machr(timeloop + 1), -Machr(timeloop + 2), 0.5 * faces.dt);
				dun_dt = secondOrderBackwardDiff(-un(timeloop), -un(timeloop + 1), -un(timeloop + 2), 0.5 * faces.dt);
				drho_dt = secondOrderBackwardDiff(-rhoCalc(timeloop), -rhoCalc(timeloop + 1), -rhoCalc(timeloop + 2), 0.5 * faces.dt);
			}
			else
			{
				dlr_dt = secondOrderBackwardDiff(Lr(timeloop), Lr(timeloop - 1), Lr(timeloop - 2), 0.5 * faces.dt);
				dMachr_dt = secondOrderBackwardDiff(Machr(timeloop), Machr(timeloop - 1), Machr(timeloop - 2), 0.5 * faces.dt);
				dun_dt = secondOrderBackwardDiff(un(timeloop), un(timeloop - 1), un(timeloop - 2), 0.5 * faces.dt);
				drho_dt = secondOrderBackwardDiff(rhoCalc(timeloop), rhoCalc(timeloop - 1), rhoCalc(timeloop - 2), 0.5 * faces.dt);
			}

			w = get_weight(faces.center({0, j}));
			p1 = dlr_dt * faces.A(j) * w / (c * src2ob.mag_r(j) * 4 * PI);
			p2 = Lr(timeloop) * faces.A(j) * w / (pow(src2ob.mag_r(j), 2) * 4 * PI);
			p3 = ((rhoCalc(timeloop) * un(timeloop) * dMachr_dt) + (1 + Machr(timeloop)) * (rhoCalc(timeloop) * dun_dt + un(timeloop) * drho_dt)) * faces.A(j) * w / (src2ob.mag_r(j) * 4 * PI);
			p4 = rhoCalc(timeloop) * un(timeloop) * ur(timeloop) * faces.A(j) * w / (pow(src2ob.mag_r(j), 2) * 4 * PI);
			pTotal(timeloop) += p1 + p2 + p3 + p4;
		}
	}
	cout << "Done.   " << endl;
}

void F1A_solver::write(size_t obs_idx)
{
	char wtf[256];
	ofstream f;
	sprintf(wtf, "%s/%s_%02d.dat", input.output_fname.c_str(), input.output_fname.c_str(), (int)obs_idx);

	f.open(wtf, std::ofstream::out);
	f << std::scientific;
	f.precision(10);
	for (size_t i = 0; i < t.get_len(); i++)
		f << setw(20) << t(i) - t(0) << setw(20) << pTotal(i) << endl;
	f.close();
}