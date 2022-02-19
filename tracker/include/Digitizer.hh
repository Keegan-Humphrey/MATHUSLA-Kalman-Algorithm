#include <iostream>
#include "globals.hh"
#include "Geometry.hh"
#include "physics.hh"
#include "par_handler.hh"

#ifndef DIGI_DEFINE
#define DIGI_DEFINE




class Digitizer{
public:
	std::vector<physics::sim_hit*> hits{};

	ParHandler* par_handler;

	int ev_num;

	int null_num = 0;

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


	



private:


}; //Digitizer




#endif
