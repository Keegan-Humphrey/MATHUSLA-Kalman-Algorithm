#include "TMinuit.h"
#include "physics.hh"
#include "LinearAlgebra.hh"
#include <iostream>
#include <fstream>

#ifndef STAT_HH
#define STAT_HH

using par = std::pair<int, int>;

class TrackFitter
{
public:
	//pars is x0, y0, z0, t0, vx, vy, vz
	static void chi2_error(int &npar, double *gin, double &f, double *par, int iflag);
	static void timeless_chi2_error(int &npar, double *gin, double &f, double *par, int iflag);
	static std::vector<physics::digi_hit *> digi_list;
	static std::vector<double> parameters;
	static std::vector<double> parameter_errors;
	const static int npar = 7;
	static double cov_matrix[npar][npar];

	bool fit(std::vector<physics::digi_hit *> _digi_list, std::vector<double> arg_guess = {})
	{
		digi_list = _digi_list;
		parameters.resize(npar);
		parameter_errors.resize(npar);

		std::vector<double> guess = arg_guess;

		TMinuit minimizer(npar);
		int ierflg = 0;
		minimizer.SetFCN(chi2_error);

		double first_step_size = 0.1;
		auto maxcalls = 500000.0;
		auto tolerance = 0.1;
		double arglist[2];
		arglist[0] = maxcalls;
		arglist[1] = tolerance;

		int quiet_mode = -1;
		int normal = 0;

		minimizer.SetPrintLevel(quiet_mode);

		double vel_y_direction = arg_guess[4]; //use the direction of the seed to help fit the track

		//we find the lowest t, and take the associated y-value and fix

		double y0_fix;

		double min_t = 9999999.99;
		for (auto hit : digi_list)
		{
			if (hit->t < min_t)
			{
				min_t = hit->t;
				y0_fix = hit->y;
			}
		}

		minimizer.mnparm(0, "x0", guess[0], first_step_size, detector::x_min, detector::x_max, ierflg);
		minimizer.mnparm(1, "y0", y0_fix, first_step_size, 0, 0, ierflg);
		minimizer.mnparm(2, "z0", guess[2], first_step_size, detector::z_min, detector::z_max, ierflg);
		minimizer.mnparm(3, "vx", guess[3], first_step_size, -1.0 * constants::c, constants::c, ierflg);
		if (vel_y_direction >= 0)
			minimizer.mnparm(4, "vy", guess[4], first_step_size, 0.0, constants::c, ierflg);
		if (vel_y_direction < 0)
			minimizer.mnparm(4, "vy", guess[4], first_step_size, -1.0 * constants::c, 0.0, ierflg);
		minimizer.mnparm(5, "vz", guess[5], first_step_size, -1.0 * constants::c, constants::c, ierflg);
		minimizer.mnparm(6, "t0", guess[6], first_step_size, 0.0, pow(10.0, 10.0), ierflg);

		minimizer.FixParameter(1);

		minimizer.mnexcm("MIGRAD", arglist, 2, ierflg);

		for (int k = 0; k < npar; k++)
		{
			minimizer.GetParameter(k, parameters[k], parameter_errors[k]);
		}

		minimizer.mnemat(&cov_matrix[0][0], npar);

		for (int k = 0; k < parameter_errors.size(); k++)
			parameter_errors[k] = TMath::Sqrt(cov_matrix[k][k]);

		//GETTTING STATUS:
		double fmin = 0.0;
		double fedm = 0.0;
		double errdef = 0.0;
		int npari = 0;
		int nparx = 0;
		/////////////////////////////////////////////////////
		int istat = 0; //this is the one we really care about
		// 0 - no covariance
		// 1 - not accurate, diagonal approximation
		// 2 - forced positive definite
		// 3 - full accurate matrix--succesful convergence

		minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);

		return (istat >= 2) ? true : false;
	}

	double rand_guess()
	{
		return 0.0;
	}

}; //class TrackFinder

class VertexFitter
{
public:
	static void nll(int &npar, double *gin, double &f, double *par, int iflag); //NEGATIVE LOG LIKLIHOOD FOR FITTING VERTICES
	static std::vector<physics::track *> track_list;
	static std::vector<double> parameters;
	static std::vector<double> parameter_errors;
	const static int npar = 4;
	static bool bad_fit;
	static double cov_matrix[npar][npar];
	static double _merit;
	//	double merit(){return _merit;}

	double merit()
	{
		double chi2 = 0.0;

		for (auto track : track_list)
		{
			chi2 += track->vertex_residual(parameters);
		}

		double ndof = static_cast<double>(npar * (track_list.size() - 1));
		ndof = ndof > 1. ? ndof : 1.0;
		return chi2 / ndof;
	}

	bool fit(std::vector<physics::track *> _track_list, std::vector<double> arg_guess = {})
	{

		track_list = _track_list;

		bad_fit = false;

		parameters.resize(npar);
		parameter_errors.resize(npar);

		std::vector<double> guess = arg_guess;

		TMinuit minimizer(npar);
		int ierflg = 0;
		minimizer.SetFCN(nll);

		double first_step_size = 0.010;
		auto maxcalls = 500000.0;
		auto tolerance = 0.01;
		double arglist[2];
		arglist[0] = maxcalls;
		arglist[1] = tolerance;

		int quiet_mode = -1;
		int normal = 0;

		minimizer.SetPrintLevel(quiet_mode);

		double min_y = detector::y_min - 10.0 * units::cm;
		double max_y = detector::y_max + 10.0 * units::cm;

		minimizer.mnparm(0, "x", guess[0], first_step_size, 0, 0, ierflg);
		minimizer.mnparm(1, "y", guess[1], first_step_size, 0, 0, ierflg);
		minimizer.mnparm(2, "z", guess[2], first_step_size, 0, 0, ierflg);
		minimizer.mnparm(3, "t", guess[3], first_step_size, 0, 0, ierflg);

		minimizer.mnexcm("MIGRAD", arglist, 2, ierflg);

		//GETTTING STATUS:
		double fmin = 0.0;
		double fedm = 0.0;
		double errdef = 0.0;
		int npari = 0;
		int nparx = 0;
		int istat = 0; //this is the one we really care about

		minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);
		_merit = fmin;

		//while (ierflg) minimizer.mnexcm("MIGRAD", arglist ,2,ierflg);
		std::ofstream file;
		file.open("print.txt", std::ios_base::app);

		file << " Parameters are " << std::endl;
		for (int k = 0; k < npar; k++)
		{
			minimizer.GetParameter(k, parameters[k], parameter_errors[k]);
			file << parameters[k] << ", ";
		}
		file << std::endl;
		//minimizer.mnemat(&cov_matrix[0][0], npar);

		return (istat >= 2) ? true : false;
	}

}; //class VertexFitter

#endif
