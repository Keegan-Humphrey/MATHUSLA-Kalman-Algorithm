#include "globals.hh"
#include "Geometry.hh"
#include "physics.hh"
#include "par_handler.hh"
#include <TRandom.h>
#include <TRandom3.h>
#include <iostream>
#include <random>

class NoiseMaker{

public:
	NoiseMaker(std::vector<physics::digi_hit*>& digis);
	static void preDigitizer(); //creates vector of all detIDs in the detector, only needs to be run once to save time
	std::vector<physics::digi_hit*> return_digis(); //returns noise digihits to digitizer to be added
	const static bool ts = false;//trouble shoot, turns on and off some messages 
	inline static bool run; //determines if the noisemaker will run or not, set as true if noise_hz in par_card > 0
private:
	void get_real_hits(std::vector<physics::digi_hit*> digis); //grabs information from digi hits about timing and location of real hits
	std::vector<double> event_timing(detID id, int hit_q); //decides the timing of the new hits for each detID
	bool get_detID_specific_hit_times(std::vector<double>* times, detID id);//gets times of real hits for specific detID	
	void make_digis(std::vector<double> times, detID id);//creates the digihits to be added
	std::vector<double> set_hit_location(detID id); //decides the location of the digihits
	
	TRandom3 hit_generator;
	std::vector<physics::digi_hit*> digi_hits; //vector of digihits created by noisemaker
	struct real_hit_times{
		detID id;
		std::vector<double> hit_times;
	}; //structure associates hit times with the correct detID
	std::vector<real_hit_times> real_hits; //vector of structures to organize real hit times and locations

	
	//static values, most used in preDigitizer()	
	static double window; //window of inserting hits after start of each event
	inline static double hits_per_second; //integer number of hits per second
	inline static double rate_of_hits;//hits_per_second divided by the correct number of nanoseconds
	inline static Geometry* _geometry;	
        inline static std::vector<detID> detID_list; //list of all detIDs in the detector
	inline static int detID_q; //total quantity of detIDs
	//static functions
        static void layer_detIDs(std::vector<detID>& _detID_list); //adds the detIDs in regular horizontal layers
	static void wall_detIDs(std::vector<detID>& _detID_list); //adds the detIDs in the vertical wall
	static void floor_detIDs(std::vector<detID>& _detID_list); //adds the detIDs in the floor 	
	

};
