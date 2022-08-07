#include<iostream>
#include <TLorentzVector.h>
#include "Geometry.hh"
#include "globals.hh"
#include "TFitter.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "LinearAlgebra.hh"
#include <Eigen/Dense>

#ifndef PHYSICS_DEFINE
#define PHYSICS_DEFINE

using Vector = vector::Vector;


namespace physics{

	//defines detector hit
	class sim_hit{
	public:

        template <typename tree_manager>
        sim_hit(tree_manager* tm, int n){
            index = n;
            x = (*(tm->sim_hit_x))[n];
            y = (*(tm->sim_hit_y))[n];
            z = (*(tm->sim_hit_z))[n];
            e = (*(tm->sim_hit_e))[n];
            t = (*(tm->sim_hit_t))[n];
            px = (*(tm->sim_hit_px))[n];
            py = (*(tm->sim_hit_py))[n];
            pz = (*(tm->sim_hit_pz))[n];
            particle_energy = (*(tm->sim_hit_particleEnergy))[n];
            pdg_id = (*(tm->sim_hit_particlePdgId))[n]; //[0];
            track_id = (*(tm->sim_hit_G4TrackId))[n]; //[0];
            particle_parent_trackid = (*(tm->sim_hit_G4ParentTrackId))[n]; //[0];
        }

		sim_hit(int index, double x, double y, double z, double t, double e){
			this->index = index;
			this->x = x;
			this->y = y;
			this->z = z;
			this->t = t;
			this->e = e;
		}

        void SetMomentum(Vector momentum){px = momentum.x; py = momentum.y; pz = momentum.z;}
        Vector GetParticleMomentum(){ return Vector(px, py, pz); }
		std::size_t index;
		double x;
		double y;
		double z;
		double t;
		double e;
        double px, py, pz; //momentum of particle which made the hit
        int track_id;
        int pdg_id;
        double particle_energy;
        double particle_parent_trackid;
		detID det_id;

	}; //sim



	class digi_hit{
	public:
		detID det_id;
		std::size_t index;
		double x, ex;
		double y, ey;
		double z, ez;
		double t;
		double e;
		double et = detector::time_resolution;
        double px, py, pz; // momentum of particle which made the hit
        double particle_mass;
        double particle_energy;
        int pdg;
        long int min_track_id = 9999999999;

        Vector PosVector(){ return Vector(x, y, z); }

		std::vector<sim_hit*> hits;

		void AddHit(sim_hit* hit){
            hits.push_back(hit);
            if (hit->track_id < min_track_id){
                min_track_id = hit->track_id;
                pdg = hit->pdg_id;
                particle_energy = hit->particle_energy;
                px = hit->px;
                py = hit->py;
                pz = hit->pz;
            }
        }



	}; //digi

	class vertex{
	public:
		double x, y, z, t;
		std::vector<int> track_indices;
	        TMatrixD _CovMatrix;
        	TMatrixD CovMatrix(){return _CovMatrix;}

		std::vector<Eigen::VectorXd> q_s;
		std::vector<Eigen::MatrixXd> D_s;
		std::vector<Eigen::VectorXd> q_f;
		std::vector<Eigen::MatrixXd> D_f;

	        double cos_opening_angle = -2.0;

        	double _merit;

	        double merit(){ return _merit; }

        	void merit(double merit){_merit = merit;}

        //set cov matrix
        template<typename matrix>
        void CovMatrix(matrix mat, int size);

		vertex(std::vector<double> pars, double cos_angle = -2){
			x = pars[0];
			y = pars[1];
			z = pars[2];
			t = pars[3];
            cos_opening_angle = cos_angle;
		}

		vertex(){}
	};

	class track{
    private:
        TMatrixD _CovMatrix;

	public:
        track(){}
		std::size_t index;
		double vx, evx;
		double vy, evy;
		double vz, evz;
		double x0, y0, z0;
		double ex0, ey0, ez0;
		double t0, et0;
	        int first_layer;
		std::vector<int> hits_to_drop = {};
		//std::vector<int> _holes;
		std::vector<int> _missing_layers;
	        std::vector<int> expected_layers;
		std::vector<int> king_move_inds;
		std::vector<double> chi_f;
        	std::vector<double> chi_s;
		std::vector<std::vector<double>> estimate_list;
		bool with_t;
		std::vector<double> x_scats;
		std::vector<double> z_scats;
		Eigen::MatrixXd P_s;
		Eigen::MatrixXd D_s;
		Eigen::MatrixXd D_f;
		Eigen::VectorXd q_s;
		Eigen::VectorXd q_f;

        //CovarianceMatrix CovMatrix(){return _CovMatrix;}
        TMatrixD CovMatrix(){return _CovMatrix;}

        //list of component hits
		std::vector<digi_hit*> hits;

        //add hit to list of hits
        //DOES NOT DO ANY RE-FITTING until used by TrackFitter
	void AddHit(physics::digi_hit* hit){ hits.push_back(hit); }
	//void set_holes(std::vector<int> layers){_holes = layers; };

        void SetExpectedLayers(std::vector<int> lyrs){expected_layers = lyrs;}

        void missing_layers(std::vector<int> layers){_missing_layers = layers; };

        //sets covariance matrix [arg mat] to more friendly form defined by the arg size
        template<typename matrix>
        void CovMatrix(matrix mat, int size);

        //constructor to automatically set the parameters and associated errors in declaration
        track(std::vector<double> params, std::vector<double> par_errors);

        //returns track parameters in vector (when called with no arguments)
        //see poorly named overload below
        std::vector<double> parameters();

        //overload to SET track parameters to arg pars
        void parameters(std::vector<double> pars);

        //SETS the paramter errros to arg epars
        void par_errors(std::vector<double> epars);

        //chi2 calculated for current parameters using all the component hits
        double chi2();

        //do I really have to comment what this is.. look at the name of the function
        double chi2_per_dof();

        //track beta
        double beta();

        //track beta error
        double beta_err();

        //2D residual (# std deviations in a given detector plane, ignores time)
        template<typename ahit>
        double untimed_residual(ahit* hit);

        //2D distance to hit in a given detector element (ignores time)
        template<typename ahit>
        double distance_to_hit(ahit* hit);

        //full 3D spatial residual at time of hit
        template<typename ahit>
        double residual(ahit* hit);

        //1D time difference using fixed y-layer
        template<typename ahit>
        double time_difference(ahit* hit);

        //1D time residual using fixed y-layer
        template<typename ahit>
        double time_residual(ahit* hit);


        //returns sorted std::vector of indices of layers in which a track has hits
        std::vector<int> layers();

         //returns the number of distinct layers in which a track has hits
        int nlayers();

        //gives 3D spatial position of a track where it would pass through a layer at the arg y
        Vector Position_at_Y(double y);

        //gives velocity vector
        Vector VelVector();

        //unit vector in direction of track
        Vector direction();

        //returns x0, y0, z0 as 3-vector
        Vector P0Vector();

        //returns track position at GLOBAL time t
        Vector position(double t);

        //gives distance between point and the track at GLOBAL time t
        double distance_to(Vector point, double t);

        //error in above distance_to calculation
        double err_distance_to(Vector point, double t);

        //gives full spatial residual for track to a vertex with arg params
        double vertex_residual(std::vector<double> params);

        //timed distance of closest approach
        double closest_approach(track* tr2);

        //closest approach of two tracks position of closest approach
        Vector closest_approach_midpoint(track* tr2);
        Eigen::VectorXd ca_midpoint_kalman(track* tr2);

        //angle the track direction forms with straight line from ip to initial point
        double cos_angle_from_ip();

        // calculate shortest dist. from track to a spatial point, by default the IP
        double shortDistance();



	}; //track

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //NOTE: TEMPLATE FUNCTIONS NEED TO BE DECALARED IN HEADER
    template<typename ahit>
    double track::untimed_residual(ahit* hit){
        double track_delta_t = (hit->y - y0)/vy;
        double track_x = x0 + vx*track_delta_t;
        double track_z = z0 + vz*track_delta_t;
        double ex2 = (hit->ex)*(hit->ex);
        double ez2 = (hit->ez)*(hit->ez);
        double res2 = (track_x - hit->x)*(track_x - hit->x)/ex2 + (track_z - hit->z)*(track_z - hit->z)/ez2;

        return TMath::Sqrt(  res2  );
    }

    template<typename ahit>
    double track::distance_to_hit(ahit* hit){
        double track_delta_t = (hit->y - y0)/vy;
        double track_x = x0 + vx*track_delta_t;
        double track_z = z0 + vz*track_delta_t;
        double res2 = (track_x - hit->x)*(track_x - hit->x) + (track_z - hit->z)*(track_z - hit->z);
        return TMath::Sqrt(  res2  );


    }

    template<typename ahit>
    double track::residual(ahit* hit){
        double track_delta_t = hit->t - t0;
        Vector track_position = Vector(x0, y0, z0) + Vector(vx, vy, vz).Scale(track_delta_t);

//	std::cout << "track residual new is " <<  track_position - hit->PosVector() << std::endl;
//	std::cout << "errors are " << Vector(hit->ex*hit->ex, hit->ey*hit->ey, hit->ez*hit->ez ) << std::endl;

        return (track_position - hit->PosVector()).Magnitude( Vector(hit->ex*hit->ex, hit->ey*hit->ey, hit->ez*hit->ez ) );
    }

    template<typename ahit>
    double track::time_difference(ahit* hit){
        
        double t_track = (hit->y - y0)/vy;
        double t_hit = hit->t - t0;
    
        return TMath::Abs(t_track - t_hit);
    }

    template<typename ahit>
    double track::time_residual(ahit* hit){
        
        double t_track = (hit->y - y0)/vy;
        double t_hit = hit->t - t0;
    
        return TMath::Abs(t_track - t_hit)/hit->et;
    }

    template<typename matrix>
    void track::CovMatrix(matrix mat, int size){

        _CovMatrix.ResizeTo(size, size);
            for (int i = 0; i < size; i++){
                for (int j = i; j < size; j++){

                    _CovMatrix[i][j] = mat[i][j];
                    _CovMatrix[j][i] = _CovMatrix[i][j];
                }
            }

    }

    template<typename matrix>
    void vertex::CovMatrix(matrix mat, int size){

        _CovMatrix.ResizeTo(size, size);
            for (int i = 0; i < size; i++){
                for (int j = i; j < size; j++){
                
                    _CovMatrix[i][j] = mat[i][j];
                    _CovMatrix[j][i] = _CovMatrix[i][j];
                }
            }   

    }



  

}; // physics



#endif
