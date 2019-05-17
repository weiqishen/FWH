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
	c_infty = sqrt(input.gamma * input.R_gas * input.T_static);//c_infty=sqrt(yRT)
	rho_infty = input.p_static / (input.R_gas * input.T_static);//rho_infty=p/RT
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
		cout << "Calculating the terms of pressure at observer location " << i + 1 << " of " << microphone.n_oberver << " as given in Farassat's 1A formulation..." << endl;

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
	cout << " -> Calculating the distance between source and microphone " << obs_id + 1 << " ... " << flush;

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
	return faces.tau(0) + src2ob.mag_r.get_max() / c_infty;
}

double F1A_solver::get_final_signal_time()
{
	return faces.tau(faces.tau.get_len() - 1) + src2ob.mag_r.get_min() / c_infty;
}

void F1A_solver::calc_pressure_term(double endt, double begint)
{
	//at source time
	ndarray<double> Lr_tau;	// Lr = Li*\hat{r}/|r| -> Li = p'\hat{n} and r = radiation vector from source faces to microphone location
	ndarray<double> rho_star_tau; //rho*=\rho_\infty+p'/c_\infty^2
	ndarray<double> Machr_tau; // M*\hat{r}/|r| of flow -> M = Mach number and r = radiation vector from source faces to mic location
	ndarray<double> un_tau;	// u*\hat{n} of flow -> Normal velocity
	ndarray<double> ur_tau;	// u*\hat{r}/|r| of flow -> velocity in radiation direction
	ndarray<double> L_tau(3), velocity_tau(3);	//Li = p\hat{n}

	//temporary variables
	size_t counter; //global time counters
	double tauCalc; //interpolate source time
	bool flag;

	//at interpolated source time
	double dlr_dt, lrCalc;
	double dMachr_dt, MachrCalc;
	double dun_dt, unCalc, urCalc;
	double drho_star_dt, rho_starCalc;
	double w, p1, p2, p3, p4;

	cout << " -> Calculating variables at interpolated source time... " << endl;

	//initialize global arrays
	//observer time array
	t.setup((size_t)((endt - begint) / faces.dt) + 1);
	for (size_t i = 0; i < t.get_len(); i++)
		t(i) = begint + i * faces.dt;

	//final pressure array
	pTotal.setup(t.get_len());
	pTotal = 0.;

	//initialize element local arrays
	Lr_tau.setup(faces.tau.get_len());
	rho_star_tau.setup(faces.tau.get_len());
	Machr_tau.setup(faces.tau.get_len());
	un_tau.setup(faces.tau.get_len());
	ur_tau.setup(faces.tau.get_len());

	//for each source
	for (size_t j = 0; j < faces.n_eles; j++)
	{
		cout << setprecision(3) << (double)j / faces.n_eles * 100 << "%   \r" << flush;

		//at source time
		for (size_t i = 0; i < faces.tau.get_len(); i++)
		{
			//p'\hat{n}
			L_tau(0) = (faces.data({j, 4, i}) - input.p_static) * faces.normal({0, j});
			L_tau(1) = (faces.data({j, 4, i}) - input.p_static) * faces.normal({1, j});
			L_tau(2) = (faces.data({j, 4, i}) - input.p_static) * faces.normal({2, j});

			//rho*=\rho_\infty+p'/c_\infty^2
			rho_star_tau(i) = rho_infty + (faces.data({j, 4, i}) - input.p_static) / pow(c_infty,2);
			//p'\hat{n}*\hat{r}/|r|, project normal stress on distance dir
			Lr_tau(i) = inner_product(L_tau.get_ptr(), L_tau.get_ptr(3), src2ob.r.get_ptr({0, j}), 0.) / src2ob.mag_r(j);

			velocity_tau(0) = faces.data({j, 1, i});
			velocity_tau(1) = faces.data({j, 2, i});
			velocity_tau(2) = faces.data({j, 3, i});

			//u_r=\hat{u}*\hat{r}/|r|
			ur_tau(i) = inner_product(velocity_tau.get_ptr(), velocity_tau.get_ptr(3), src2ob.r.get_ptr({0, j}), 0.) / src2ob.mag_r(j);
			//\hat{u}*\hat{n} normal velocity
			un_tau(i) = inner_product(velocity_tau.get_ptr(), velocity_tau.get_ptr(3), faces.normal.get_ptr({0, j}), 0.);
			//mach_r=u_r/c_infty
			Machr_tau(i) = ur_tau(i) / c_infty;
		}

		counter = 0;
		//at interpolated source time
		for (size_t timeloop = 0; timeloop < t.get_len(); timeloop++) //for each observer time step
		{
			tauCalc = t(timeloop) - src2ob.mag_r(j) / c_infty; //the interpolated src time

			//calculate interpolated variables
			flag = false;
			while (counter < faces.tau.get_len() - 1)
			{
				if (tauCalc >= faces.tau(counter) && tauCalc < faces.tau(counter + 1))
				{
					if (tauCalc != faces.tau(counter)) //use 2 neighbouring soure time
					{
						dlr_dt = calc_slope(Lr_tau(counter), Lr_tau(counter + 1), faces.dt);
						dMachr_dt = calc_slope(Machr_tau(counter), Machr_tau(counter + 1), faces.dt);
						dun_dt = calc_slope(un_tau(counter), un_tau(counter + 1), faces.dt);
						drho_star_dt = calc_slope(rho_star_tau(counter), rho_star_tau(counter + 1), faces.dt);
					}
					else
					{
						if (counter > 0)
						{
							dlr_dt = calc_slope(Lr_tau(counter - 1), Lr_tau(counter + 1), 2 * faces.dt);
							dMachr_dt = calc_slope(Machr_tau(counter - 1), Machr_tau(counter + 1), 2 * faces.dt);
							dun_dt = calc_slope(un_tau(counter - 1), un_tau(counter + 1), 2 * faces.dt);
							drho_star_dt = calc_slope(rho_star_tau(counter - 1), rho_star_tau(counter + 1), 2 * faces.dt);
						}
						else
						{
							dlr_dt = calc_slope(Lr_tau(counter), Lr_tau(counter + 1), faces.dt);
							dMachr_dt = calc_slope(Machr_tau(counter), Machr_tau(counter + 1), faces.dt);
							dun_dt = calc_slope(un_tau(counter), un_tau(counter + 1), faces.dt);
							drho_star_dt = calc_slope(rho_star_tau(counter), rho_star_tau(counter + 1), faces.dt);
						}
					}

					lrCalc = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, Lr_tau(counter), Lr_tau(counter + 1));
					MachrCalc = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, Machr_tau(counter), Machr_tau(counter + 1));
					unCalc = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, un_tau(counter), un_tau(counter + 1));
					urCalc = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, ur_tau(counter), ur_tau(counter + 1));
					rho_starCalc = linearInterpolate(faces.tau(counter), faces.tau(counter + 1), tauCalc, rho_star_tau(counter), rho_star_tau(counter + 1));

					w = get_weight(faces.center({0, j}));
					p1 = dlr_dt / (c_infty * src2ob.mag_r(j));
					p2 = lrCalc / pow(src2ob.mag_r(j), 2);
					p3 = ((rho_starCalc * unCalc * dMachr_dt) + (1 + MachrCalc) * (rho_starCalc * dun_dt + unCalc * drho_star_dt)) / src2ob.mag_r(j);
					p4 = rho_starCalc * unCalc * urCalc / pow(src2ob.mag_r(j), 2);
					pTotal(timeloop) += (p1 + p2 + p3 + p4) * (faces.A(j) * w) / (4. * PI);

					flag = true;
					break;
				}
				counter++;
			}
			if (flag == false)
				Fatal_Error("cant find interpolated source time");
		}
	}
	cout << "Done.   " << endl;
}

double F1A_solver::get_weight(double x)
{
	double eps = 1e-10;
	if (input.endcap_avg)
	{
		if (x < (input.endcap_x(0) - eps))
			return 1.0;
		for (size_t i = 0; i < input.n_endcaps - 1; i++)
		{
			if (fabs(x - input.endcap_x(i)) <= eps)
				return 1.0 / (double)input.n_endcaps;
			else if (x > input.endcap_x(i) + eps && x < input.endcap_x(i + 1) - eps)
				return (double)(input.n_endcaps - 1 - i) / (double)input.n_endcaps;
		}
		if (fabs(x - input.endcap_x(input.n_endcaps - 1)) <= eps)
			return 1.0 / (double)input.n_endcaps;

		return 0.0;
	}
	else
	{
		return 1.0;
	}
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