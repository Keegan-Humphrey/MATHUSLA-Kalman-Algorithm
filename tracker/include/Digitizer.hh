#include <iostream>
#include "globals.hh"
#include "Geometry.hh"
#include "physics.hh"
#include "par_handler.hh"
#include <TRandom3.h>

#ifndef DIGI_DEFINE
#define DIGI_DEFINE




class Digitizer{
public:
	std::vector<physics::sim_hit*> hits{};

	ParHandler* par_handler;

	int ev_num;

	long long int seed;

	int dropped_hits = 0;
	int floor_wall_hits = 0;

	int null_num = 0;

	TRandom3 generator;
	TRandom3 drop_generator;

	std::vector<physics::digi_hit*> Digitize();
	Geometry* _geometry;

	void AddHit(physics::sim_hit *hit){
		hit->det_id = (_geometry->GetDetID(hit));
		if (!hit->det_id.IsNull()) {
			hits.push_back(hit);
		}
		else null_num++;

	}

	void clear(){
		
		for (auto hit : hits) delete hit;
		hits.clear();
	
	}

	Digitizer(){ _geometry = new Geometry; }

	~Digitizer() 
	{
		delete _geometry;
		for (auto p : hits) delete p;
	}


	
	void InitGenerators(){
		TRandom3 _generator;

        	int seed_size = static_cast<int>(1e9);

	        seed = par_handler->par_map["seed"];
        	seed = seed == -1 ? rand()*rand()*rand() % seed_size : seed;

	        if (par_handler->par_map["debug"] == 1) {
        	        std::cout << "Digi seed is: " << seed << std::endl;
        	}

	        _generator.SetSeed(seed);

	        // Now we throw out hits in the floor and wall to simulate reduced detector efficiency
		std::vector<physics::digi_hit*> digis_not_dropped;

	        //TRandom drop_generator;
        	TRandom3 _drop_generator;
	        _drop_generator.SetSeed( rand()*rand()*rand() % rand() );

		generator = _generator;
		drop_generator = _drop_generator;
	}


private:


}; //Digitizer




#endif
