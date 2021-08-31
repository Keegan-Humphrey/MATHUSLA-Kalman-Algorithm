#include<cstdlib>
#include<vector>
#include <fstream>
#include <string>

#ifndef UNITS_HH
#define UNITS_HH

namespace units{
	const double mm = 10.0;
	const double cm = 1.00;
	const double ns = 1.00;
	const double MeV = 1.00;
};

namespace detector{
	using namespace units;
	const double ip_x = 0.0;
	const double ip_y = 0.0;
	const double ip_z = 0.0;

	//specifies the bottom and top y position of every layer
	const std::vector<std::vector<double>> LAYERS_Y={{6003.0*cm, 6006.0*cm},  //layer 0
 												 	{6106.0*cm, 6109.0*cm}, //layer 1
													{6209.0*cm, 6212.0*cm}, //layer 2
 													{8003.0*cm, 8006.0*cm}, //layer 3
 													{8106.0*cm, 8109.0*cm}, //layer 4
 													{8503.0*cm, 8506.0*cm}, //layer 5
 													{8606.0*cm, 8609.0*cm}, //layer 6
 													{8709.0*cm, 8712.0*cm}, //layer 7
 													{8812.0*cm, 8815.0*cm}, //layer 8
													{8915.0*cm, 8918.0*cm} };  //layer 9
	const int n_layers = 10;

	const std::vector<std::vector<double>> MODULE_X = { 	{-4950.0*cm, -4050.0*cm},
													{-3950.0*cm, -3050.0*cm},
													{-2950.0*cm, -2050.0*cm},
													{-1950.0*cm, -1050.0*cm},
													{-950.0*cm, -50.0*cm},
													{50.0*cm, 950.0*cm},
													{1050.0*cm, 1950.0*cm},
													{2050.0*cm, 2950.0*cm},
													{3050.0*cm,  3950.0*cm},
													{4050.0*cm,  4950.0*cm} };

	const std::vector<std::vector<double>> MODULE_Z = {{7000.0*cm, 7900.0*cm},
													{8000.0*cm, 8900.0*cm},
													{9000.0*cm, 9900.0*cm},
													{10000.0*cm, 10900.0*cm},
													{11000.0*cm, 11900*cm},
													{12000.0*cm, 12900.0*cm},
													{13000.0*cm, 13900.0*cm},
													{14000.0*cm, 14900.0*cm},
													{15000.0*cm, 15900.0*cm},
													{16000.0*cm, 16900.0*cm} };

	const int n_modules = 100;
	const double scintillator_length = 450.0*units::cm;
	const double scintillator_width = 4.50*units::cm;
	const double scintillator_height = 3.0*units::cm;
	const double scintillator_thickness = 2.0*units::cm; //Just the sensitive part
	const double time_resolution = 1.0*units::ns;

	const double x_min = MODULE_X[0][0];
	const double y_min = LAYERS_Y[0][0];
	const double z_min = MODULE_Z[0][0];

	const double x_max = MODULE_X[MODULE_X.size()-1][1];
	const double y_max = LAYERS_Y[LAYERS_Y.size()-1][1];
	const double z_max = MODULE_Z[MODULE_Z.size()-1][1];

	//FLOOR TILE WIDTHS

	const double floor_x_width = 50.0*units::cm;
	const double floor_z_width = 50.0*units::cm;

    //WALL TILE WIDTHS

    const double wall_x_width = 50.0*units::cm;
    const double wall_y_width = 50.0*units::cm;

    //WALL PARAMETERS

    const double wall_gap = 1.0*units::cm; //gap on each side of wall
    const double wall_height = 2000.0*units::cm;
    const double wall_start_y = 6000.0*units::cm; //min y value of wall (casing included)

    //FOR statistics.hh ONLY - NEW MIN Z WITH WALL
    const double z_min_wall = z_min - wall_gap - scintillator_height;

};

namespace constants{

	const double c = 29.97*units::cm/units::ns;
	const double optic_fiber_n = 1.500; //estimate for the optical fiber index of refraction

};

namespace cuts{

	//digi cuts and constants
	const double digi_spacing = 20.0*units::ns;
	const double SiPM_energy_threshold = 0.65*units::MeV;

	//seeding
	const double seed_ds2 = 5.0; //sigma
	const double seed_residual = 10.0; //sigma

	//tracking
	const double residual_drop = 12.0; //sigma
	const double residual_add = 12.0; //sigma
	const double track_chi2 = 5.0;
	const int track_nlayers = 3; //4; // used to be 3
	const int nseed_hits = 4;
	const double time_difference_drop = 12.0; //sigma
	const double seed_time_difference = 10.0; //ns
	const int ntrack_hits = 4;
	const double distance_to_hit = 75.0*units::cm;

	//kalman tracker
        const double kalman_chi_s = 150.0;
        const double kalman_chi_add = 200.0;
        const double kalman_track_chi = 15.0;
        const std::vector<double> kalman_v_add = {0.800000,1.200000};
        const std::vector<double> kalman_v_drop = {0.900000,1.100000};

	//merging
        const double merge_cos_theta = 0.998;
        const double merge_distance = 25.0*units::cm;

	//cleaning step
	const double chi2_add = 10.0;
	const double chi2_drop = 16.0;
	const int cleaning_nhits = 6;

	//vertexing
        const double seed_closest_approach = 100.0*units::cm;
        const double vertex_chi2 = 15.0;
        const double closest_approach_add = 100.0*units::cm;
        const double kalman_vertex_chi_add = 100000.0;
        const double kalman_vertex_chi = 100.0;

	//run options
        const int start_ev = 0;
        const int end_ev = 50;

    //digi hit cuts for floors and wall
    const std::vector<bool> include_floor = { true, true, true }; //ith index for ith floor from bottom

    const bool include_wall = true;
    const double wall_y_cut = detector::wall_start_y + detector::wall_height; //all digi hits above this will be thrown out
};


namespace kalman{

	const double sigma_ms_p = 5.01; // [rad MeV]
        const double p = 500.0; // [MeV] representative momentum

//        const std::vector<double> v_cut_add = {-2.0,2.0};
//        const std::vector<double> v_cut_drop = {-2.0,2.0};
        const std::vector<double> v_cut_add = {0.0,2.0};
        const std::vector<double> v_cut_drop = {0.5,1.5};

	const double pull_cut_add = 100.0;
	const double pull_cut_drop = 100.0;
}



#endif
