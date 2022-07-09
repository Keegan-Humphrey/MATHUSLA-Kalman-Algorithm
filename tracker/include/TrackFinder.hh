#include "physics.hh"
#include "globals.hh"
#include "LinearAlgebra.hh"
#include <TMath.h>
#include <TRandom.h>
#include "Geometry.hh"
#include <iostream>
#include <fstream>
#include "par_handler.hh"

#ifndef TF_DEFINE
#define TF_DEFINE

//#include "kalman.hh"
//#include <Eigen/Dense>

class seed
{
public:
	double score;

	std::pair<physics::digi_hit *, physics::digi_hit *> hits;
	seed() {}
	seed(physics::digi_hit *hit1, physics::digi_hit *hit2) : hits(hit1, hit2)
	{
	}

	vector::Vector VelVector()
	{
		auto P1 = hits.first->PosVector();
		auto P2 = hits.second->PosVector();
		double dt = hits.second->t - hits.first->t;

		return (P2 - P1).Scale(1.0 / dt);
	}

	std::vector<double> guess()
	{
		auto P1 = hits.first->PosVector();
		auto velocity = VelVector();

		return {P1.x, P1.y, P1.z, velocity.x, velocity.y, velocity.z, hits.first->t};
	}

	template <typename hit>
	double time_difference(hit AHit)
	{

		auto P1 = hits.first->PosVector();
		auto vel = VelVector();
		double seed_t = (AHit->y - P1.y) / vel.y;

		double hit_t = AHit->t - hits.first->t;

		return TMath::Abs(seed_t - hit_t) / AHit->et;
		;
	}

	template <typename hit>
	double timeless_residual(hit AHit)
	{

		//calculate residual from the track using the two hits
		//(assuming a straight line between them)

		auto P1 = hits.first->PosVector();
		auto velocity = VelVector();

		//NOW we use the TIMELESS residual!!!
		//this gives us a measure of how far the hit is from the track when the track would be in that layer

		auto HitP = AHit->PosVector();
		double delta_t = (HitP.y - P1.y) / velocity.y;
		auto expectedPos = P1 + velocity.Scale(delta_t);
		auto errMetric = vector::Vector(pow(AHit->ex, 2), pow(AHit->ey, 2), pow(AHit->ez, 2));

		return (expectedPos - HitP).Magnitude(errMetric);
	}

	template <typename hit>
	double distance_to_hit(hit AHit)
	{

		auto P1 = hits.first->PosVector();
		auto velocity = VelVector();

		//NOW we use the TIMELESS residual!!!
		//this gives us a measure of how far the hit is from the track when the track would be in that layer

		auto HitP = AHit->PosVector();
		double delta_t = (HitP.y - P1.y) / velocity.y;
		auto expectedPos = P1 + velocity.Scale(delta_t);

		return (expectedPos - HitP).Magnitude();
	}
};

//#include "kalman-test.hh" // must be after seed class definition

class TrackFinder
{
public:
	std::ofstream file;

	std::vector<physics::digi_hit *> hits;
	std::vector<physics::digi_hit *> hits_k;
	std::vector<seed> seeds;
	std::vector<seed> seeds_k;
	std::vector<physics::track *> tracks;
	std::vector<physics::track *> tracks_k;
	std::vector<physics::track *> tracks_k_m;
	std::vector<double> local_chi_f;
	std::vector<double> local_chi_s;
	std::vector<int> failure_reason;

	std::vector<physics::digi_hit *> track_pts;
	std::vector<physics::digi_hit *> unused_hits;
	std::vector<physics::digi_hit *> good_hits;
	std::vector<physics::digi_hit *> undropped_hits;

	void CalculateMissingHits(Geometry *geo);
	void CalculateHoles(Geometry *geo);
	void clear()
	{
		seeds.clear();
		for (auto hit : hits)
			delete hit;
		hits.clear();
		for (auto tr : tracks)
			delete tr;
		tracks.clear();

		seeds_k.clear();
//		for (auto hit : hits_k)
//			delete hit;
		// already taken care of by hits delete statement

		hits_k.clear();

		for (auto tr : tracks_k)
			delete tr;
		tracks_k.clear();
		// for (auto tr : tracks_k_m)
		// 	delete tr;
		tracks_k_m.clear();

		local_chi_f.clear();
		local_chi_s.clear();
	}

	void clear_vecs()
	{

		track_pts.clear();
		unused_hits.clear();
		undropped_hits.clear();
		good_hits.clear();
	}

	seed first_seed;
	//	kalman_do kft;
	ParHandler* par_handler;

	int first_n_to_delete = 0;
	void Seed();
	void FindTracks();
	void FindTracks_kalman();
	void FindTracks_all();
	void CleanTracks();
	void Reseed(bool);
	void CheckSeeds();
	void MergeTracks();
	void MergeTracks_k();

	//	bool seed_unused(seed current_seed);
	//	bool extract_good_hits();

	bool seed_unused(seed current_seed)
	{
		//std::ofstream file;
		//file.open("print.txt", std::ios_base::app);

		bool seed_hit_in = false;
		for (auto hit : hits_k)
		{
			if (hit->index == current_seed.hits.first->index)
			{
				seed_hit_in = true;
				break;
			}
		}
		if (!seed_hit_in)
		{
			file << "seed's first hit was already used by another track!!!!" << '\n';
		}
		//file.close();
		return seed_hit_in;
	}

	double _score();
	double c_score(physics::digi_hit *hit1, physics::digi_hit *hit2)
	{

		double dx = hit1->x - hit2->x;
		double dy = hit1->y - hit2->y;
		double dz = hit1->z - hit2->z;
		double dt = hit1->t - hit2->t;

		return TMath::Abs((dx * dx + dy * dy + dz * dz) / (constants::c * constants::c) - dt * dt);
	} //c_score

	double c_score(const seed &s)
	{
		auto val = c_score(s.hits.first, s.hits.second);
		;

		return val;
	}

	void score_seeds()
	{
		//std::ofstream file;
		//file.open("print.txt", std::ios_base::app);

		for (auto seed : seeds)
			seed.score = c_score(seed);

		int i;
		for (auto seed : seeds_k)
		{
			auto P1 = seed.hits.first->PosVector();
			auto P2 = seed.hits.second->PosVector();

			double dr = (P2 - P1).Magnitude() / constants::c; // [ns]

			seeds_k[i].score = c_score(seed) / (dr * dr); // these will now be ordered by relative scores

			i++;
		}
		//file.close();
	}

	int min_seed()
	{ //sorts by c_compatability score

		int min_index = -1;
		double min_val = cuts::seed_ds2;
		int j = 0;
		for (auto seed : seeds)
		{
			if (seed.score < min_val)
			{
				min_index = j;
				min_val = seed.score;
			}

			j++;
		}

		return min_index;
	}

	int min_seed_k()
	{ //sorts by c_compatability score
		//std::ofstream file;
		//file.open("print.txt", std::ios_base::app);

		int min_index = -1;
		double min_val = 1e6; // large value to be overwritten by first seed
		int j = 0;
		for (auto seed : seeds_k)
		{
			if (seed.score < min_val)
			{

				min_index = j;
				min_val = seed.score;
			}

			j++;
		}

		//file.close();

		return min_index;
	}
}; //TrackFinder

#endif
