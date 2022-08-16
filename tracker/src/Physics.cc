#include<iostream>
#include <TLorentzVector.h>
#include "Geometry.hh"
#include "globals.hh"
#include "TFitter.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "LinearAlgebra.hh"
#include "physics.hh"
#include <Eigen/Dense>


using Vector = vector::Vector;

namespace physics {


	track::track(std::vector<double> params, std::vector<double> par_errors){
			x0 = params[0]; ex0 = par_errors[0];
			y0 = params[1]; ey0 = par_errors[1];
			z0 = params[2]; ez0 = par_errors[2];
			vx = params[3]; evx = par_errors[3];
			vy = params[4]; evy = par_errors[4];
			vz = params[5]; evz = par_errors[5];
			if (params.size() == 7) {
				t0 = params[6];
				et0 = par_errors[6];
				with_t = true;
			}
			calculateAngularValues();
			calculateAngularError();

	}
	
	track::track(std::vector<double> params, std::vector<double> par_errors,bool angular){
		if(angular){
			x0 = params[0]; ex0 = par_errors[0];
			y0 = params[1]; ey0 = par_errors[1];
			z0 = params[2]; ez0 = par_errors[2];
			theta = params[3]; etheta = par_errors[3];
			phi = params[4]; ephi = par_errors[4];
			mbeta = params[5]; ebeta = par_errors[5];
			if (params.size() == 7) {
				t0 = params[6];
				et0 = par_errors[6];
				with_t = true;
			}
			calculateCartesianValues();
			calculateCartesianErrors();
		}
	
	}

	void track::calculateAngularValues(){
		theta = cmath::atan2(cmath::sqrt(vx*vx+vx*vz),vy);
		phi = cmath::atan2(vz,vx);
		mbeta = beta();
		}

	void track::calculateAngularError(){
		etheta = 0;
		ephi = 0;
		ebeta = beta_err();
		//TODO: calculate proper errors for these values
	}

	void track::calculateCartesianValues(){
		vx = mbeta*constants::c*cmath::sin(theta)*cmath::cos(phi);
		vy = mbeta*constants::c*cmath::cos(theta);
		vz = mbeta*constants::c*cmath::sin(theta)*cmath::sin(phi);
	}
	
	void track::calculateCartesianError(){
		evx = 0;
		evy = 0;
		evz = 0;
		//TODO: calculate proper errors for these values
	}

	std::vector<double> track::parameters(){
			std::vector<double> p = { x0, y0, z0, vx, vy, vz };
			return p;
	}
	
	std::vector<double> track::angularparameters(){
			std::vector<double> p = { x0, y0, z0, theta, phi, mbeta }
			return p;
	}

	void track::parameters(std::vector<double> pars){
			x0 = pars[0];
			y0 = pars[1];
			z0 = pars[2];
			vx = pars[3];
			vy = pars[4];
			vz = pars[5];
			if (pars.size() == 7) {
				t0 = pars[6];
				with_t = true;
			}
			calculateAngularValues();
		}

	void track::angularparameters(std::vector<double> pars){
			x0 = pars[0];
			y0 = pars[1];
			z0 = pars[2];
			theta = pars[3];
			phi = pars[4];
			mbeta = pars[5];
			if (pars.size() == 7) {
				t0 = pars[6];
				with_t = true;
			}
			calculateCartesianValues();
		}

	void track::par_errors(std::vector<double> epars){
			ex0 = epars[0];
			ey0 = epars[1];
			ez0 = epars[2];
			evx = epars[3];
			evy = epars[4];
			evz = epars[5];
			if (epars.size() == 7) {
				et0 = epars[6];
				with_t = true;
			}
			calculateAngularError();
		}
	
	void track::angular_par_errors(std::vector<double> epars){
			ex0 = epars[0];
			ey0 = epars[1];
			ez0 = epars[2];
			etheta = epars[3];
			ephi = epars[4];
			ebeta = epars[5];
			if (epars.size() == 7) {
				et0 = epars[6];
				with_t = true;
			}
			calculateCartesianError();
		}
	

	double track::chi2(){
			double chi_2 = 0.0;
			double t;
			double res_t, res_x, res_z;
			for (auto hit : hits){
				t = (hit->y - y0)/vy;

				res_x = ((x0 + vx*t) - hit->x)/hit->ex;
				res_z = ((z0 + vz*t) - hit->z)/hit->ez;

				if  (!with_t) chi_2 +=  res_x*res_x + res_z*res_z;
				else {
					res_t = ((t0 + t) - hit->t)/hit->et;
					chi_2 +=  res_x*res_x + res_z*res_z + res_t*res_t;

					//std::cout << "zero pos are " << x0 << ", " << t0 << ", " << z0 << ", " << y0 << std::endl;
					//std::cout << "zero res are " << res_x << ", " << res_t << ", " << res_z << std::endl;

				}
			}

			return chi_2;
		}

    double track::chi2_per_dof(){
      	int n_track_params = 6;
      	int ndof = (4.0*hits.size() - n_track_params );
      	if (ndof < 1) ndof = 1;
      	return chi2()/ndof; //FIX ME
      }

    double track::beta(){
      	return TMath::Sqrt( vx*vx + vy*vy + vz*vz  )/constants::c;
      }

    double track::beta_err(){
        double norm = beta()*(constants::c * constants::c);
        std::vector<double> derivatives = { 0., 0., vx/norm, vy/norm, vz/norm, 0.};

      	double error = 0.0;

        for (int i = 0; i < derivatives.size(); i++){
            for (int j = 0; j < derivatives.size(); j++){

                error += derivatives[i]*derivatives[j]*_CovMatrix[i][j];

            }
        }

        return TMath::Sqrt(error);
      }


    std::vector<int> track::layers(){
    	std::vector<int> layer_indices;
    	for (auto hit : hits){
    		bool new_layer = true;
    		int layer_index = hit->det_id.layerIndex;
    		for (int layer_n : layer_indices){
    			if (layer_n == layer_index){
    				new_layer = false;
    				break;
    			}
    		}
	    	if (new_layer) layer_indices.push_back(layer_index);
    	}
        std::sort(layer_indices.begin(), layer_indices.end());
    	return layer_indices;
    } //nlayers


    int track::nlayers(){
    	//returns the number of layers that a track has hits in
    	return layers().size();
    } //nlayers


    Vector track::Position_at_Y(double y){
    	double delta_t = (y-y0)/vy;
    	return Vector(x0, y0, z0) + Vector(vx, vy, vz).Scale(delta_t);
    }

    Vector track::VelVector(){
    	return Vector(vx, vy, vz);
    }

    Vector track::direction(){
    	double velocity = beta()*constants::c;
    	return VelVector().Scale(1.0/velocity);
    }

    Vector track::P0Vector(){
    	return Vector(x0, y0, z0);
    }

    Vector track::position(double t){ //global time t
    	double dt = t-t0;
    	return P0Vector() + VelVector().Scale(dt);
    }

    double track::distance_to(Vector point, double t){
    	auto pos = position(t);
    	return (pos - point).Magnitude();
    }

    double track::err_distance_to(Vector point, double t){

    	std::vector<double> derivatives = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //partial derivates of position_at
    	//note, we don't include the derivative with respect to y. y0 is fixed, and so it isnt included in the covariance matrix
    	auto x = point.x;
    	auto y = point.y;
    	auto z = point.z;

    	double dist = distance_to(point, t);

    	derivatives[0] = ((t-t0)*vx - x + x0)/dist;
    	//derivatives[1] = detector::scintillator_height;
    	derivatives[1] = ((t-t0)*vz - z + z0)/dist;
    	derivatives[2] = (t-t0)*((t-t0)*vx - x + x0)/dist;
    	derivatives[3] = (t-t0)*((t-t0)*vy - y + y0)/dist;
    	derivatives[4] = (t-t0)*((t-t0)*vz - z + z0)/dist;
    	derivatives[5] = -1.0*(vx*((t-t0)*vx - x + x0) + vy*((t-t0)*vy - y + y0) + vz*((t-t0)*vz - z + z0))/dist;

    	//now we calculate the actual error
    	double error = 0.0;

    	for (int i = 0; i < derivatives.size(); i++){
    		for (int j = 0; j < derivatives.size(); j++){

    			error += derivatives[i]*derivatives[j]*_CovMatrix[i][j];

    		}
    	}

    	return TMath::Sqrt(error);
    }


    double track::vertex_residual(std::vector<double> params){

    	double x = params[0];
    	double y = params[1];
    	double z = params[2];
    	double t = params[3];

    	return distance_to(Vector(x, y, z), t)/err_distance_to(Vector(x, y, z), t);
    }


    //SEE docs/Closest_Approach_Of_Parametric_Vectors.pdf for the math
    double track::closest_approach(track* tr2){
    	using namespace vector;
        if (tr2->index == index) return 0.00;

        Vector rel_v = tr2->VelVector() - this->VelVector();
    	double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

    	Vector displacement = this->P0Vector() - tr2->P0Vector();

    	double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

    	Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
    	Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);

    	return (pos1 - pos2).Magnitude();
    }


    //SEE docs/Closest_Approach_Of_Parametric_Vectors.pdf for the math
    Vector track::closest_approach_midpoint(track* tr2){

    	using namespace vector;
        if (tr2->index == index) return P0Vector();

        Vector rel_v = tr2->VelVector() - this->VelVector();
    	double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

    	Vector displacement = this->P0Vector() - tr2->P0Vector();

    	double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

    	Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
    	Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);
    
    	return (pos1 + pos2).Scale(0.5);
    }

    Eigen::VectorXd track::ca_midpoint_kalman(track* tr2){

        using namespace vector;
//        if (tr2->index == index) return P0Vector();

        Vector rel_v = tr2->VelVector() - this->VelVector();
        double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

        Vector displacement = this->P0Vector() - tr2->P0Vector();

        double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

        Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
        Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);

        Vector mid = (pos1 + pos2).Scale(0.5);

	Eigen::VectorXd mid_t(4);

	mid_t[0] = mid.x;
	mid_t[1] = mid.y;
	mid_t[2] = mid.z;
	mid_t[3] = t_ca;

	return mid_t;
    }


    double track::cos_angle_from_ip(){

        using namespace detector;

        Vector ip_direction = Vector(x0 - ip_x, y0 - ip_y, z0 - ip_z );
        double track_vel = beta()*constants::c;

        return (ip_direction ^ VelVector())/(track_vel*ip_direction.Magnitude());
    }

    // calculate shortest dist. from point to line
    double track::shortDistance()
    {
        //      std::cout << "222: " << line_point2 << " 111: " << line_point1 << "\n";
        Vector point = Vector(detector::ip_x, detector::ip_y, detector::ip_z );
        vector::Vector AB = vector::Vector(vx, vy, vz);
        vector::Vector AC = point - vector::Vector(x0, y0, z0);
        double area = vector::Vector(AB * AC).Magnitude();
        double CD = area / AB.Magnitude();
        return CD;
    }
    





}; //namespcae physics
