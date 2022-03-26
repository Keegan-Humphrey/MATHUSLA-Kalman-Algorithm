#include "Digitizer.hh"
#include "physics.hh"
#include "globals.hh"
#include <TRandom.h>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <array>
#include <iostream>
#include <vector>
#include <sys/time.h>
//#include "units.hh"
//DIGITIZATION ALGORITHM
namespace physics{
	bool time_sort(physics::sim_hit *hit1, physics::sim_hit *hit2) {return (hit1->t < hit2->t);}
};

std::vector<physics::digi_hit*> Digitizer::Digitize(){

	//this is the vector of digi_hits we will return at the end of the function
	std::vector<physics::digi_hit*> digis;


	//looping through each detector ID
	std::vector<physics::sim_hit*> current_hits;
	std::vector<physics::sim_hit*> current_remaining_hits = hits;
	std::vector<physics::sim_hit*> next_remaining_hits;

	while (current_remaining_hits.size() > 0){

		//current detector id which we are working in
		auto current_id = (current_remaining_hits[0])->det_id;


		//taking out all hits with the same detector id to be digitized, leaving the remaing for the next iteration
		for (auto hit : current_remaining_hits){
			if (hit->det_id.IsNull()) continue;
			if (hit->det_id == current_id){ current_hits.push_back(hit);}
			else {next_remaining_hits.push_back(hit);}

		}

		//time sorting current hits

		std::sort(current_hits.begin(), current_hits.end(), &physics::time_sort);


		// going through all hits until they are either all added to digis, or dropped

		while (current_hits.size() > 0){

			std::vector<physics::sim_hit*> used_hits;
			std::vector<physics::sim_hit*> unused_hits;

			double t0 = (current_hits[0])->t;

			double e_sum = 0;

			for (auto hit : current_hits){
				if ( hit->t < t0 + cuts::digi_spacing ){
					e_sum += hit->e;
					used_hits.push_back(hit);
				} else { unused_hits.push_back(hit);}
			}
            
            // ignoring all hits in ignored floors/walls or above wall y cut
            
            auto current_center = _geometry->GetCenter(current_id);
            if (current_id.isFloorElement){
                if (cuts::include_floor[current_id.layerIndex] != true){
                    current_hits.erase(current_hits.begin());
                    continue;
                }
            }
            if (current_id.isWallElement){
                if (cuts::include_wall != true || current_center[1] > cuts::wall_y_cut){
                    current_hits.erase(current_hits.begin());
                    continue;
                }
            }
            
			if (e_sum > cuts::SiPM_energy_threshold){
				physics::digi_hit* current_digi = new physics::digi_hit();
				current_digi->det_id = current_id;
				for (auto hit : used_hits){current_digi->AddHit(hit);}
				current_digi->index = ( digis.size() );
				digis.push_back(current_digi);
				current_hits = unused_hits;
			}  else { current_hits.erase(current_hits.begin());}

		}

		//resetting all the sorting vectors, and assigning the next remianing hits to the next iteration for current remaining
		current_remaining_hits.clear();
		current_remaining_hits = next_remaining_hits;
		next_remaining_hits.clear();
		current_hits.clear();

	}

	//at this point, all of the digi_hits in the digi_vector have the hits which will make them up. However, they don't have any of their energy, position, or timing information added.
	//Below, we compute the energy, time, and position of all of the digi hits
	//We incoorporate the time and position smearing into this calculation as well

	int counter = 0;

	struct timeval curTime;
	gettimeofday(&curTime, NULL);
	long int micro_sec = curTime.tv_usec;
	srand( micro_sec );
//	srand( time(NULL) + ev_num ); // incase called more than once per second across events

	TRandom generator;
	//TRandom3 generator;

//	generator.SetSeed( rand()*rand()*rand() % rand() );
//	int seed = 14557409;
//	int seed = 897765236; // bad eff seed
//	int seed = -182399494; // good eff seed
//	int seed = rand()*rand()*rand() % rand();

        int seed = par_handler->par_map["seed"];
        seed = seed == -1 ? rand()*rand()*rand() % rand() : seed;

        if (par_handler->par_map["debug"] == 1) {
                std::cout << "Digi seed is: " << seed << std::endl;
	}

	generator.SetSeed(seed);

	// Now we throw out hits in the floor and wall to simulate reduced detector efficiency
	std::vector<physics::digi_hit*> digis_not_dropped;

	//std::cout << par_handler->par_map["scint_efficiency"];
	//std::cout << 1.0 / par_handler->par_map["scint_efficiency"] << std::endl;

	TRandom drop_generator;
	drop_generator.SetSeed( rand()*rand()*rand() % rand() );

	for (auto digi : digis) {
		// Rndm() uniformly samples (0,1) so 1 time in every scint_efficiency samples we don't satisfy this condition
		if (drop_generator.Rndm() > 1.0 / par_handler->par_map["scint_efficiency"]) {
			digis_not_dropped.push_back(digi);
		}
		else std::cout << "dropped a hit" << std::endl;
	}

	digis.clear();
	digis = digis_not_dropped;
	digis_not_dropped.clear();

	// now manage hits in the floor and wall
	for (auto digi : digis){
	//	std::cout << counter++ << std::endl;
		auto current_id = digi->det_id;

		auto center = _geometry->GetCenter(current_id);
		auto layer = _geometry->layer_list[0];
		auto long_direction_index = layer->long_direction_index;
		auto uncertainty = layer->uncertainty();

		if (current_id.isFloorElement){
			//if (current_id.isFloorElement) std::cout << digi->hits.size() << std::endl;
			uncertainty = _geometry->_floor.uncertainty(current_id.layerIndex);
        } else if (current_id.isWallElement){
            uncertainty = _geometry->_wall.uncertainty();
		} else {
			layer = _geometry->layer_list[current_id.layerIndex];
			long_direction_index = layer->long_direction_index;
			uncertainty = layer->uncertainty();
		}

		double e_sum = 0;
		double long_direction_sum = 0.0;
		double t_sum = 0;

		for (auto hit : digi->hits){
			e_sum += hit->e;
			t_sum += hit->t * hit->e;

			if (long_direction_index == 0){
				long_direction_sum += hit->x * hit->e;
			} else {
				long_direction_sum += hit->z * hit->e;
			}
		}

		//std::cout << "getting energy weighting" << std::endl;

	//	std::cout << "uncertainty:" << std::endl;
	//	for (int count = 0; count < uncertainty.size(); count++){
	//		std::cout << "count: " << count << "  element: " << uncertainty[count] << std::endl;
	//	}
	//	std::cout << "center:" << std::endl;
	//	for (int count = 0; count < center.size(); count++){
	//		std::cout << "count: " << count << "  element: " << center[count] << std::endl;
	//	}
		digi->e = e_sum;
		digi->t = t_sum/e_sum;
		digi->y = center[1];
		digi->ey = uncertainty[1];
		digi->ex = uncertainty[0];
		digi->ez = uncertainty[2];

	//	std::cout << "done" << std::endl;

		//note: et is the same for all of them and is set in the digi class defintion
		if (current_id.isFloorElement || current_id.isWallElement){
			digi->x = center[0];
			digi->z = center[2];
        } else {
	    	if (long_direction_index == 0){
				digi->x = long_direction_sum/e_sum;
				digi->z = center[2];
			} else {
				digi->z= long_direction_sum/e_sum;
				digi->x = center[0];
			}
		}

		//TIME AND POSITION SMEARING!!!!!!!!!!!!!!!
		//we see the random number generator with a number that should be completly random:
		//the clock time times the layer index times the number of digis


		digi->t += generator.Gaus(0.0, digi->et);

		if (current_id.isFloorElement || current_id.isWallElement) {

			continue;
		}

		if (long_direction_index == 0) {
			double smeared_x = digi->x + generator.Gaus(0.0, digi->ex);
			if (!(_geometry->GetDetID(smeared_x, digi->y, digi->z) == current_id) ){
				if (smeared_x > center[0]) { smeared_x = center[0] + (layer->widths())[0]/2.0 - 0.5*units::cm; }
				else {smeared_x = center[0] - (layer->widths())[0]/2.0 + 0.5*units::cm; }
			}

			digi->x = smeared_x;

		} else if (long_direction_index == 1){
			double smeared_z = digi->z + generator.Gaus(0.0, digi->ez);
			if (!(_geometry->GetDetID(digi->x, digi->y, smeared_z) == current_id) ){
				if (smeared_z > center[2]) {smeared_z = center[2] + (layer->widths())[1]/2.0 - 0.5*units::cm;}
				else {smeared_z = center[2] - (layer->widths())[1]/2.0 + 0.5*units::cm; }
			}
			digi->z = smeared_z;
		}

		if ( !(_geometry->GetDetID(digi->x, digi->y, digi->z) == current_id) ){
			std::cout << "Warning!!! Smearing function error--digi was smeared to be outside of known detector element!!" << std::endl;
		}

	}

	//setting digi indices
	int k = 0;
	for (auto digi : digis) {
		digi->index = k++;
	}

	return digis;



}
