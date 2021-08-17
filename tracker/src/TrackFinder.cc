#include "TrackFinder.hh"
#include "globals.hh"
#include "physics.hh"
#include <TMath.h>
#include <TRandom.h>
#include "TMinuit.h"
#include "statistics.hh"
#include "Geometry.hh"
#include "kalman.hh"
#include "kalman-test.hh"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

void TrackFinder::Seed()
{
//	std::ofstream file;
//	file.open("print.txt", std::ios_base::app);

	// if (file.fail()) std::cout << "open didn't work 1" << std::endl;
	if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

	seeds.clear();
	seeds = {};

	for (int first = 0; first < hits.size(); first++)
	{
		for (int second = first + 1; second < hits.size(); second++)
		{

			int layer1 = ((hits[first])->det_id).layerIndex;
			int layer2 = ((hits[second])->det_id).layerIndex;

			if (layer1 == layer2)
				continue;

			double ds = c_score(hits[first], hits[second]);

			if (ds > cuts::seed_ds2)
				continue;

			seeds.push_back(seed(hits[first], hits[second]));
			seeds_k.push_back(seed(hits[first], hits[second]));

			file << "seed:" << hits[first]->index << "," << hits[second]->index << '\n';
		} //"second" loop
	}	  //"first" loop

	score_seeds();

	//for (auto seed : seeds) file << c_score(seed.hits.first, seed.hits.second) << std::endl;

	file << "Seeding finished!" << '\n';

	file.close();
} //TF::Seed

void TrackFinder::FindTracks()
{
//	std::ofstream file;
	file.open("print.txt", std::ios_base::app);

	if (file.fail()) std::cout << "open didn't work 2" << std::endl;
	if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

	file << "----------------------------------------------" << std::endl;
	file << "--------------Begin old tracking--------------" << '\n';
	file << "----------------------------------------------" << std::endl;

	file.close();

	if (seeds.size() == 0)
		return; //no seeds found in initial seeding, will be retried with <c travel

	//we take the first seed now

	int index = 0;
	int total_hits = hits.size();
	bool iterate = true;
	int j = 0;
	int MAX_ITS = 25;

	while (iterate)
	{

		if (seeds.size() == 0)
			return;
		if (hits.size() == 0)
			return;

	        file.open("print.txt", std::ios_base::app);

		int min_index = min_seed();
		auto current_seed = seeds[min_index];

		//std::ofstream file;

		if (file.fail()) std::cout << "open didn't work 3" << std::endl;
		if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

		file << "------------------------------" << std::endl;
		file << "Current_seed:" << current_seed.hits.first->index << "," << current_seed.hits.second->index << '\n';
		file << "------------------------------" << std::endl;

		seeds.erase(seeds.begin() + min_index); //delete the seed so that it isn't used again

		std::vector<physics::digi_hit *> track_pts;
		std::vector<physics::digi_hit *> unused_hits;

		//now we use the seed::residual function to get the residual from the seed of all of the available hits.
		//if the residual is good, we add it to the track_pts vector

		for (auto hit : hits)
		{
			// file << " TR " << current_seed.timeless_residual(hit) << " DH " << current_seed.distance_to_hit(hit) << " TD " << current_seed.time_difference(hit) << std::endl;
			if ((current_seed.timeless_residual(hit) < cuts::seed_residual or current_seed.distance_to_hit(hit) < cuts::distance_to_hit) and current_seed.time_difference(hit) < cuts::seed_time_difference)
			{
				track_pts.push_back(hit);
				// file << "this hit-> seed:"<< hit->index << '\n';
			}
			else
			{

				unused_hits.push_back(hit);
				// file << "this hit-> unused_hits:"<< hit->index << '\n';
			}
		}

		//we check that there are at least cuts::nseed_hits hits in the track_pts
		//if not, we move on to the next seed

		if (track_pts.size() < cuts::nseed_hits)
		{
			file << "seed failed!" << '\n';
			file.close();
			continue;
		}

		// file << "now first fit!"<< '\n';

		//at this point, all we have done is erase the seed, and none of the hits have been modified
		TrackFitter fitter;
		auto track_status = fitter.fit(track_pts, current_seed.guess());

		if (track_status == false) {
			file.close();
			continue; //fit failed
		}

		file << "Fit parameters are " << fitter.parameters[3] << " , " << fitter.parameters[4] << " , " << fitter.parameters[5] << " , " << std::endl;
		auto current_track = new physics::track(fitter.parameters, fitter.parameter_errors);
		for (auto hit : track_pts)
			current_track->AddHit(hit);

		//first pass track is now cconstructed

		std::vector<physics::digi_hit *> second_unused_hits;

		for (auto hit : unused_hits)
		{
			// file << " TR " << current_track->untimed_residual(hit) << " TD " << current_track->time_difference(hit) << std::endl;
			if ((current_track->untimed_residual(hit) < cuts::residual_add or current_track->distance_to_hit(hit) < cuts::distance_to_hit) and current_track->time_difference(hit) < cuts::seed_time_difference)
			{
				current_track->AddHit(hit);
				// file << "this hit also -> 1st track " << hit->index << '\n';
			}
			else
			{
				second_unused_hits.push_back(hit);
				// file << "this hit -> 2nd_unused " << hit->index << '\n';
			}
		}

		//NOW WE REFIT AGAIN

		fitter.fit(current_track->hits, current_track->parameters());
		current_track->parameters(fitter.parameters);
		current_track->par_errors(fitter.parameter_errors);

		//file << " y is " << std::endl;
		//for (auto hit : current_track->hits) file << hit->y << ", ";
		//file << std::endl;

		for (int i = 0; i < fitter.parameters.size() - 1; i++)
			file << "old track parameters:" << fitter.parameters[i] << '\n';
		for (int i = 0; i < fitter.parameter_errors.size() - 1; i++)
			file << "old track parameters errors:" << fitter.parameter_errors[i] << '\n';

		// file << "now second fit!"<< '\n';

		std::vector<physics::digi_hit *> good_hits;
		for (auto hit : current_track->hits)
		{

			// file << " TR " << current_track->untimed_residual(hit) <<" TD " << current_track->time_difference(hit) << std::endl;
			if (current_track->untimed_residual(hit) > cuts::residual_drop or current_track->time_difference(hit) > cuts::time_difference_drop)
			{
				second_unused_hits.push_back(hit);
				// file << "this hit -> 2nd_unused " << hit->index << '\n';
			}
			else
			{
				good_hits.push_back(hit);
				//file << "this hit -> made the track " << hit->index << '\n';
				//file << "this hit has residual " << current_track->residual(hit) << std::endl;
			}
		}
		//		file << " chi squared is: " << current_track->chi2_per_dof() << std::endl;

		if (good_hits.size() < cuts::nseed_hits)
		{
			file << "track failed!" << good_hits.size() << '\n';
			file.close();
			continue;
		}

		current_track->hits = good_hits;

		fitter.fit(current_track->hits, current_track->parameters());
		current_track->parameters(fitter.parameters);
		current_track->par_errors(fitter.parameter_errors);
		current_track->CovMatrix(fitter.cov_matrix, fitter.npar);

		file << " chi squared is: " << current_track->chi2_per_dof() << std::endl;
		file << " n layers is: " << current_track->nlayers() << std::endl;
		file << " n hits is: " << current_track->hits.size() << std::endl;

		if (current_track->nlayers() >= cuts::track_nlayers and current_track->chi2_per_dof() < cuts::track_chi2 and current_track->hits.size() > cuts::ntrack_hits)
		{
			tracks.push_back(current_track);
			file << "track made it!" << '\n';

		}
		else
		{
			delete current_track;
			file << "track failed" << '\n';
			file.close();
			continue;
		}

		hits = second_unused_hits;

		if (seeds.size() == 0)
			iterate = false;
		if (hits.size() < cuts::nseed_hits)
			iterate = false;

		file.close();
	} //while iterate

	//MergeTracks();
	//CleanTracks();

	int total_tracked_points = 0;

	for (auto track : tracks)
		total_tracked_points += track->hits.size();

	if ((total_hits - hits.size() - total_tracked_points) != 0)
	{
		file << "total: " << total_hits << std::endl;
		file << "used: " << total_tracked_points << std::endl;
		file << "unused:" << hits.size() << std::endl;
	}

	if (j++ > MAX_ITS)
		iterate = false;

	file.close();
}

void TrackFinder::MergeTracks()
{
//	std::ofstream file;
	file.open("print.txt", std::ios_base::app);

	//at this point, all of the points have been fit to tracks. At this point, we will perform the track merging step.

	if (tracks.size() == 0 or tracks.size() == 1)
		return;

	std::vector<int> deleted_tracks = {};

	for (int first_track = 0; first_track < tracks.size(); first_track++)
	{
		for (int second_track = first_track + 1; second_track < tracks.size(); second_track++)
		{

			auto tr1 = tracks[first_track];
			auto tr2 = tracks[second_track];

			auto first_layer = tr1->layers()[0] < tr2->layers()[0] ? tr1->layers()[0] : tr2->layers()[0];

			auto d1 = tr1->direction();
			auto d2 = tr2->direction();

			int tr1_missing_hits = 0;
			int tr2_missing_hits = 0;

			auto cos_theta = d1 ^ d2;

			double distance = tr1->closest_approach(tr2);

			if (distance > cuts::merge_distance or cos_theta < cuts::merge_cos_theta)
				continue;

			//at this point, we need to check if they have a certain number of missing hits

			//if (tr1->_missing_layers.size() < 3 or tr2->_missing_layers.size() < 3) continue;
			std::vector<int> joint_missing_hit_layers = {};
			for (int i = 0; i < tr1->_missing_layers.size(); i++)
			{
				auto layer = tr1->_missing_layers[i];
				if (layer < first_layer)
					continue;
				tr1_missing_hits++;

				bool missing = true;
				for (int j = 0; j < tr2->_missing_layers.size(); j++)
				{

					auto layer2 = tr2->_missing_layers[j];

					if (layer2 == layer)
					{
						missing = false;
						continue;
					}
				}

				if (missing)
					joint_missing_hit_layers.push_back(layer);
			}

			for (int i = 0; i < tr2->_missing_layers.size(); i++)
			{
				auto layer = tr2->_missing_layers[i];
				if (layer < first_layer)
					continue;
				tr2_missing_hits++;
				bool missing = true;
				for (int j = 0; j < tr1->_missing_layers.size(); j++)
				{

					auto layer2 = tr2->_missing_layers[j];

					if (layer2 == layer)
					{
						missing = false;
						continue;
					}
				}

				for (int j = 0; j < joint_missing_hit_layers.size(); j++)
				{

					auto layer2 = joint_missing_hit_layers[j];

					if (layer2 == layer)
					{
						missing = false;
						continue;
					}
				}

				if (missing)
					joint_missing_hit_layers.push_back(layer);
			}

			bool merge = false;

			if (joint_missing_hit_layers.size() < 3 and (tr1_missing_hits > 2 or tr2_missing_hits > 2))
				merge = true;

			if (!merge)
				continue;

			for (auto hit : tr2->hits)
				tr1->AddHit(hit);

			TrackFitter fitter;

			fitter.fit(tr1->hits, tr1->parameters());
			tr1->parameters(fitter.parameters);
			tr1->par_errors(fitter.parameter_errors);

			deleted_tracks.push_back(second_track);

		} //second track
	}	  //first track

	std::vector<physics::track *> good_tracks;

	for (int k = 0; k < tracks.size(); k++)
	{

		bool add = true;
		for (int del_index : deleted_tracks)
		{
			if (del_index == k)
				add = false;
		}

		if (add)
		{
			good_tracks.push_back(tracks[k]);
		}
	}

	tracks = good_tracks;

	//at this point, the list of tracks is finalized. Now they will be indexed:
	int trackn = 0;
	for (auto tr : tracks)
	{
		tr->index = trackn++;
	}
	file.close();
}

void TrackFinder::MergeTracks_k()
{
//	std::ofstream file;
	// file.open("print.txt", std::ios_base::app);

	std::ofstream file;
	file.open("cout.txt", std::ios_base::app);

	//at this point, all of the points have been fit to tracks. At this point, we will perform the track merging step.

	if (tracks_k_m.size() == 0 or tracks_k_m.size() == 1)
		return;

	std::vector<int> deleted_tracks = {};

	file << "-------------------- Commence Merge Fit -------------" << std::endl;

	for (int first_track = 0; first_track < tracks_k_m.size(); first_track++)
	{

		for (int second_track = first_track + 1; second_track < tracks_k_m.size(); second_track++)
		{

			file << "--------------- tracks " << first_track << ", " << second_track << std::endl;

			auto tr1 = tracks_k_m[first_track];
			auto tr2 = tracks_k_m[second_track];

			// std::cout << "tr1_hole:" << tr1->_holes.size() << '\n';
			// std::cout << "tr2_hole:" << tr2->_holes.size() << '\n';

			bool mergebool = false;


			auto first_layer = tr1->layers()[0] < tr2->layers()[0] ? tr1->layers()[0] : tr2->layers()[0];

			auto d1 = tr1->direction();
			auto d2 = tr2->direction();

			int tr1_missing_hits = 0;
			int tr2_missing_hits = 0;

			auto cos_theta = d1 ^ d2;

			double distance = tr1->closest_approach(tr2);

			file << "opening angle is " << cos_theta << std::endl;
			file << "closest approach is " << distance << std::endl;

			if (distance < cuts::merge_distance and cos_theta > cuts::merge_cos_theta){
				mergebool = true;
			}

			if ((tr1->_holes.size() >= 3 or tr2->_holes.size() >= 3) and (distance < 150 and cos_theta > 0.95)){
				mergebool = true;
			}

			if ((tr1->_holes.size() >= 2 or tr2->_holes.size() >= 2) and (distance < 100 and cos_theta > 0.99)){
				mergebool = true;
			}
			//at this point, we need to check if they have a certain number of missing hits
			file << "tr1->_missing_layers: "<< tr1->_missing_layers.size() << '\n';
			file << "tr2->_missing_layers: "<< tr2->_missing_layers.size() << '\n';

			//if (tr1->_missing_layers.size() < 3 or tr2->_missing_layers.size() < 3) continue;
			std::vector<int> joint_missing_hit_layers = {};
			for (int i = 0; i < tr1->_missing_layers.size(); i++)
			{
				auto layer = tr1->_missing_layers[i];
				if (layer < first_layer)
					continue;
				tr1_missing_hits++;

				bool missing = true;
				for (int j = 0; j < tr2->_missing_layers.size(); j++)
				{

					auto layer2 = tr2->_missing_layers[j];

					if (layer2 == layer)
					{
						missing = false;
						continue;
					}
				}

				if (missing)
					joint_missing_hit_layers.push_back(layer);
			}

			for (int i = 0; i < tr2->_missing_layers.size(); i++)
			{
				auto layer = tr2->_missing_layers[i];
				if (layer < first_layer)
					continue;
				tr2_missing_hits++;
				bool missing = true;
				for (int j = 0; j < tr1->_missing_layers.size(); j++)
				{

					auto layer2 = tr2->_missing_layers[j];

					if (layer2 == layer)
					{
						missing = false;
						continue;
					}
				}

				for (int j = 0; j < joint_missing_hit_layers.size(); j++)
				{

					auto layer2 = joint_missing_hit_layers[j];

					if (layer2 == layer)
					{
						missing = false;
						continue;
					}
				}

				if (missing)
					joint_missing_hit_layers.push_back(layer);
			}

			bool merge = false;
			// std::cout << "joint_missing_hit_layers.size()" << joint_missing_hit_layers.size() << '\n';
			// std::cout << "tr1_missing_hits: "<< tr1_missing_hits  << '\n';
			// std::cout << "tr2_missing_hits: "<< tr2_missing_hits  << '\n';

			if (joint_missing_hit_layers.size() < 3 and (tr1_missing_hits > 2 or tr2_missing_hits > 2))
				merge = true;

			//if (!merge)
			//	continue;
			if (!mergebool){continue;}

			seed artificial_seed = seed(tr1->hits[0], tr1->hits[1]);

			for (auto hit : tr2->hits)
				tr1->AddHit(hit);

			kalman_track kft;

			kft.finding = false;
			kft.dropping = false;
			kft.kalman_all(tr1->hits, &artificial_seed);

			deleted_tracks.push_back(second_track);

			if (kft.status != 2) // if fit fails, delete second track and don't change first
				continue;

			tr1->hits = kft.added_hits;

			tr1->x_scats = kft.x_scat;
			tr1->z_scats = kft.z_scat;

			tr1->chi_f = kft.chi_f;
			tr1->chi_s = kft.chi_s;

			tr1->estimate_list = kft.x_s_list;
			tr1->P_s = kft.P_s0;

			for (auto ind : tr2->king_move_inds) tr1->king_move_inds.push_back(ind);

		} //second track
	}	  //first track

	std::vector<physics::track *> good_tracks;

	for (int k = 0; k < tracks_k_m.size(); k++)
	{

		bool add = true;
		for (int del_index : deleted_tracks)
		{
			if (del_index == k)
				add = false;
		}

		if (add)
		{
			good_tracks.push_back(tracks_k_m[k]);
		}
	}

	// std::cout << "tracks before merge: " << tracks_k_m.size() << std::endl;

	tracks_k_m = good_tracks;

	// std::cout << "tracks after merge: " << tracks_k_m.size() << std::endl;

	//at this point, the list of tracks is finalized. Now they will be indexed:
	for (int trackn=0; trackn<tracks_k_m.size(); trackn++)
	{
		tracks_k_m[trackn]->index = trackn;
	}
	file.close();
}

void TrackFinder::CleanTracks()
{
//	std::ofstream file;
	file.open("print.txt", std::ios_base::app);

	//we clean the tracks with the following critera:
	//-

	if (tracks.size() < 2)
		return;

	for (int first_track = 0; first_track < tracks.size(); first_track++)
	{
		for (int second_track = first_track + 1; second_track < tracks.size(); second_track++)
		{
			auto tr1 = tracks[first_track];
			auto tr2 = tracks[second_track];

			if (tr1->hits.size() < cuts::cleaning_nhits or tr2->hits.size() < cuts::cleaning_nhits)
				continue;

			file << "clean me" << std::endl;

		} //second track
	}	  //first track
	file.close();
}

std::vector<physics::digi_hit *> TrackFitter::digi_list = {};
std::vector<double> TrackFitter::parameters = {};
std::vector<double> TrackFitter::parameter_errors = {};
double TrackFitter::cov_matrix[TrackFitter::npar][TrackFitter::npar];

void TrackFitter::chi2_error(int &npar, double *gin, double &f, double *pars, int iflag)
{

	double x0 = pars[0];
	double y0 = pars[1];
	double z0 = pars[2];
	double vx = pars[3];
	double vy = pars[4];
	double vz = pars[5];
	double t0 = pars[6];

	double error = 0.0;

	for (auto hit : TrackFitter::digi_list)
	{
		double t = (hit->t - t0);
		double dt = (hit->y - y0) / vy;
		double expected_x = x0 + dt * vx;
		double expected_z = z0 + dt * vz;
		double _ex = (expected_x - hit->x) / hit->ex;
		double _ez = (expected_z - hit->z) / hit->ez;
		double _et = (t - dt) / hit->et;
		error += (_ex * _ex + _ez * _ez) + _et * _et;
	}

	f = error;
}

void TrackFinder::CalculateHoles(Geometry *geo){
	for (auto track : tracks_k_m)
	{
		// std::cout << "--new track:" << '\n';
		//file << "new track " << std::endl;
		std::vector<int> layers = track->layers();
		std::vector<int> expected_layers;


		int layer_n = track->layers()[0];
		for (auto layer_lims : detector::LAYERS_Y)
		{

			double y_center = (layer_lims[0] + layer_lims[1]) / 2.0;
			auto track_position = track->Position_at_Y(y_center);

			if (track_position.x > detector::x_min and track_position.x < detector::x_max)
			{
				if (track_position.z > detector::z_min and track_position.z < detector::z_max)
				{
					if (!(geo->GetDetID(track_position).IsNull()))
						expected_layers.push_back(layer_n);
				}
			}

			layer_n++;
			if (layer_n == track->layers().back()){break;}
		}

		// track->SetExpectedLayers(expected_layers);

		std::vector<int> missing_layers;

		for (auto expected_index : expected_layers)
		{

			bool missing = true; //flag to indicate if the "expected_index" for the layer is missing or not

			for (auto existing_index : layers)
			{
				if (expected_index == existing_index)
					missing = false;
			}

			bool already_counted = false;
			if (missing)
			{
				for (auto _index : missing_layers)
				{
					if (expected_index == _index)
						already_counted = true;
				}

				if (!already_counted)
				{
					//	file << expected_index << std::endl;
					missing_layers.push_back(expected_index);
				}
			} //if missing
		}
		// for (auto missedlayer : missing_layers) {
		// 	std::cout << "missed:"<< missedlayer << '\n';
		// }
		// std::cout << "***missednum: "<< missing_layers.size() << '\n';
		track->set_holes(missing_layers);
	}
}

void TrackFinder::CalculateMissingHits(Geometry *geo)
{

	for (auto track : tracks)
	{
		// std::cout << "--new track:" << '\n';
		//file << "new track " << std::endl;
		std::vector<int> layers = track->layers();
		std::vector<int> expected_layers;

		int layer_n = 0;
		for (auto layer_lims : detector::LAYERS_Y)
		{

			double y_center = (layer_lims[0] + layer_lims[1]) / 2.0;
			auto track_position = track->Position_at_Y(y_center);

			if (track_position.x > detector::x_min and track_position.x < detector::x_max)
			{
				if (track_position.z > detector::z_min and track_position.z < detector::z_max)
				{
					if (!(geo->GetDetID(track_position).IsNull()))
						expected_layers.push_back(layer_n);
				}
			}

			layer_n++;
		}

		track->SetExpectedLayers(expected_layers);

		std::vector<int> missing_layers;

		for (auto expected_index : expected_layers)
		{

			bool missing = true; //flag to indicate if the "expected_index" for the layer is missing or not

			for (auto existing_index : layers)
			{
				if (expected_index == existing_index)
					missing = false;
			}

			bool already_counted = false;
			if (missing)
			{
				for (auto _index : missing_layers)
				{
					if (expected_index == _index)
						already_counted = true;
				}

				if (!already_counted)
				{
					//	file << expected_index << std::endl;
					missing_layers.push_back(expected_index);
				}
			} //if missing
		}
		for (auto missedlayer : missing_layers) {
			// std::cout << "missed:"<< missedlayer << '\n';
		}
		track->missing_layers(missing_layers);
	}
	//calculate for kalman tracks

	for (auto track : tracks_k_m)
	{
		//file << "new track " << std::endl;
		std::vector<int> layers = track->layers();
		std::vector<int> expected_layers;

		int layer_n = 0;
		for (auto layer_lims : detector::LAYERS_Y)
		{

			double y_center = (layer_lims[0] + layer_lims[1]) / 2.0;
			auto track_position = track->Position_at_Y(y_center);

			if (track_position.x > detector::x_min and track_position.x < detector::x_max)
			{
				if (track_position.z > detector::z_min and track_position.z < detector::z_max)
				{
					if (!(geo->GetDetID(track_position).IsNull()))
						expected_layers.push_back(layer_n);
				}
			}

			layer_n++;
		}

		track->SetExpectedLayers(expected_layers);

		std::vector<int> missing_layers;

		for (auto expected_index : expected_layers)
		{

			bool missing = true; //flag to indicate if the "expected_index" for the layer is missing or not

			for (auto existing_index : layers)
			{
				if (expected_index == existing_index)
					missing = false;
			}

			bool already_counted = false;
			if (missing)
			{
				for (auto _index : missing_layers)
				{
					if (expected_index == _index)
						already_counted = true;
				}

				if (!already_counted)
				{
					//	file << expected_index << std::endl;
					missing_layers.push_back(expected_index);
				}
			} //if missing
		}

		track->missing_layers(missing_layers);
	}
}

void TrackFinder::FindTracks_kalman()
{
//	std::ofstream file;
	file.open("print.txt", std::ios_base::app);

	if (file.fail()) std::cout << "open didn't work 4" << std::endl;
	if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

	file << "----------------------------------------------" << std::endl;
	file << "--------------Begin new tracking--------------" << '\n';
	file << "----------------------------------------------" << std::endl;
	file.close();

	if (seeds_k.size() == 0)
		return; //no seeds found in initial seeding, will be retried with <c travel

	int index = 0;
	int total_hits = hits_k.size();
	bool iterate = true;
	int j = 0;
	int MAX_ITS = 25;

	while (iterate)
	{
		if (seeds_k.size() == 0)
			//return;
			break;
		if (hits_k.size() == 0)
			//return;
			break;

		int min_index = min_seed_k();
		auto current_seed = seeds_k[min_index];

		//std::ofstream file;
	        file.open("print.txt", std::ios_base::app);

		if (file.fail()) std::cout << "open didn't work 5" << std::endl;
		if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

		file << "CURRENT SEED SCORE IS: " << current_seed.score << std::endl;
		file << "------------------------------" << std::endl;
		file << "Current_seed:" << current_seed.hits.first->index << "," << current_seed.hits.second->index << '\n';
		file << "------------------------------" << std::endl;
		file.close();

		seeds_k.erase(seeds_k.begin() + min_index); //delete the seed so that it isn't used again

		bool used = !seed_unused(current_seed); // double negative! used = not unused

//		std::cout << "hits size 1: " << hits_k.size() << std::endl;

		clear_vecs();

		std::vector<int> king_move_inds;

//		std::cout << "hits size 2: " << hits_k.size() << std::endl;

		// Construct the first filter (to find good hits for seed)
		kalman_track kf_find;
		kf_find.finding = true;
		kf_find.dropping = true;
		kf_find.seed_was_used = used;
		kf_find.kalman_all(hits_k, &current_seed);

		king_move_inds = kf_find.king_move_inds;

		if (kf_find.status == -1)
		{
			failure_reason[5] += 1;
			continue;
		}

		if (kf_find.status == -2)
                {
                        failure_reason[6] += 1;
                        continue;
                }

//		good_hits = kf_find.found_hits;
		undropped_hits = kf_find.found_hits;
		unused_hits = kf_find.unadded_hits;

		int drops = -1;
		int i = 0;

		bool failed = false;

		//std::vector<kalman_track *> kfts;
//		kalman_track kft_;

//		for (int i=0; i < 3; i++)
		while (drops != 0)
		{
			file << "Drop Iteration " << i << std::endl;

			kalman_track kft_;
        	        kft_.finding = false;
                	kft_.dropping = true;
	                kft_.seed_was_used = used;

			//kft.unadded_hits = kf_find.unadded_hits;
			kft_.unadded_hits = unused_hits;

                	//kft.kalman_all(kf_find.found_hits, &current_seed);
                	//kft.kalman_all(good_hits, &current_seed);
                	kft_.kalman_all(undropped_hits, &current_seed);

			//kft_ = kft; // new
			//kfts.push_back(&kft);

			file.open("print.txt", std::ios_base::app);

			if (file.fail()) std::cout << "open didn't work 6" << std::endl;
			if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

			if (kft_.status == -1)
        	        {
                	        failure_reason[5] += 1;
				file.close();
				failed = true;
				break;
	                        //continue;
        	        }

			good_hits = kft_.added_hits;

			undropped_hits.clear(); // new

			if (kft_.status == -2)
			{
				failure_reason[6] += 1;
				file.close();
				failed = true;
				break;
//				continue;
			}

			if (good_hits.size() < cuts::track_nlayers)
			{
				failure_reason[4] += 1;
			}

			unused_hits = kft_.unadded_hits;
			unused_hits.insert(unused_hits.end(), kf_find.unadded_hits.begin(), kf_find.unadded_hits.end());

	//		std::cout << "first unused length " << unused_hits.size() << std::endl;
	//		std::cout << "second unused length " << unused_hits.size() << std::endl;

	//		file << "hits size 3: " << hits_k.size() << std::endl;

			drops = 0;

			//dropping hits
			for (int n = 0; n < good_hits.size(); n++)
			{
				Eigen::VectorXd v(3);
				v << kft_.v_s_list[n][0], kft_.v_s_list[n][1], kft_.v_s_list[n][2];

				file << "beta is " << v.norm() / constants::c << std::endl;
				file << "beta drop " << !(cuts::kalman_v_drop[0] < v.norm() / constants::c
							&& v.norm() / constants::c < cuts::kalman_v_drop[1]) << std::endl;
				if (kft_.chi_s[n] > cuts::kalman_chi_s
			           || !(cuts::kalman_v_drop[0] < v.norm() / constants::c && v.norm() / constants::c < cuts::kalman_v_drop[1]))
				{
					file << "dropped with chi: " << kft_.chi_s[n] << '\n';
					unused_hits.push_back(good_hits[n]);
					drops++;
				}
				else
				{
					file << "added with chi: " << kft_.chi_s[n] << '\n';
					undropped_hits.push_back(good_hits[n]);
				}
			}

			file << "number of hits dropped " << drops << std::endl;
//			std::cout << "number of hits dropped " << drops << std::endl;

			if (undropped_hits.size() < cuts::track_nlayers)
			{
				file << "track failed! with " << undropped_hits.size() << " hits" << '\n';
				failure_reason[0] += 1;
				file.close();
				failed = true;
				break;
//				continue;
			}
			file.close();

			//king_move_inds = kft_.king_move_inds;
//			for (auto ind : kft_.king_move_inds) king_move_inds.push_back(ind);

			//// end where kft used to be

			if (i == 15) {std::cout << "track i break " << std::endl; break;}

//			std::cout << "Iteration " << i << std::endl;
			i++;

			if (drops != 0) continue; // new
//			drops = 0;

		// don't do an extra fit?
//		kalman_track kft_2;
//		kft_2.finding = false;
//		kft_2.dropping = false;
//		kft_2.kalman_all(undropped_hits, &current_seed);

//		kft_ = kft_2
		//kalman_track kft_ = *kfts.back();

//		for (auto ind : kft_.king_move_inds) king_move_inds.push_back(ind);

		//if (kft_.king_move_inds.size() != 0) std::cout << "\n trackfinder first";
		//for (auto ind : kft_.king_move_inds) std::cout << ind << ", ";

		file.open("print.txt", std::ios_base::app);

		if (file.fail()) std::cout << "open didn't work 7" << std::endl;
		if (file.bad()) std::cout << "bad: " << file.bad() << std::endl;

		file << " status is " << kft_.status << std::endl;

		if (kft_.status == 0)
		{
			failure_reason[1] += 1;
			file.close();
//			continue;
			failed = true;
			break;
		}
		if (kft_.status == -1)
		{
			failure_reason[8] += 1;
			file.close();
//			continue;
			failed = true;
			break;
		}
		if (kft_.status == -2)
		{
			failure_reason[7] += 1;
			file.close();
//			continue;
			failed = true;
			break;
		}

//		file << "test 0" << std::endl;

		auto current_track = new physics::track(kft_.x_s, kft_.P_s);
		current_track->hits = undropped_hits;

//		file << "test 1" << std::endl;

		current_track->x_scats = kft_.x_scat;
		current_track->z_scats = kft_.z_scat;

//		file << "test 2" << std::endl;

		current_track->chi_f = kft_.chi_f;
		current_track->chi_s = kft_.chi_s;

//		file << "test 3" << std::endl;

		current_track->estimate_list = kft_.x_s_list;
		current_track->P_s = kft_.P_s0;

//		file << "test 4" << std::endl;

		current_track->king_move_inds = kft_.king_move_inds;

		//if (king_move_inds.size() != 0) std::cout << "\n trackfinder second";
		//for (auto ind : king_move_inds) std::cout << ind << ", ";

		static double cov_matrix[7][7];
		for (int i = 0; i < 7; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				if (i == j)
				{
					cov_matrix[i][j] = std::pow(kft_.P_s[i],2);
				}
				else
				{
					cov_matrix[i][j] = 0;
				}
			}
		}

		current_track->CovMatrix(cov_matrix, 7);

		for (auto chi : kft_.chi_f)
			local_chi_f.push_back(chi);
		for (auto chi : kft_.chi_s)
			local_chi_s.push_back(chi);

		file << "chi f is ";
		for (auto chi : kft_.chi_f)
			file << chi << " , ";
		file << std::endl;

		file << "chi s is ";
		for (auto chi : kft_.chi_s)
			file << chi << " , ";
		file << std::endl;

		double chi_sum = 0;
		for (auto chi : kft_.chi_s)
			chi_sum += chi;
		chi_sum = chi_sum / (4.0 * kft_.chi_s.size() - 6.0);

		file << "chi s sum is " << chi_sum << std::endl;

		file << " chi squared is: " << current_track->chi2_per_dof() << std::endl;
		file << " n layers is: " << current_track->nlayers() << std::endl;
		file << " n hits is: " << current_track->hits.size() << std::endl;

		if (current_track->nlayers() >= cuts::track_nlayers && chi_sum < cuts::kalman_track_chi)
		{
			tracks_k.push_back(current_track);
			file << "track made it!" << '\n';
			// std::cout << "kalman track layers: "<< '\n';
			// for (auto layer : current_track->layers()){std::cout << "layer: "<< layer << '\n';}
		}
		else
		{
			if (current_track->nlayers() < cuts::track_nlayers)
				failure_reason[2] += 1;
			if (current_track->hits.size() < cuts::ntrack_hits)
				failure_reason[3] += 1;
			file << "track failed" << '\n';
			delete current_track;
			file.close();
			failed = true;
			break;
//			continue;
		}


//		std::cout << "third unused length " << unused_hits.size() << std::endl;

		} // drop for loop

		if (failed) continue;

		hits_k = unused_hits;

		file << "seeds k size is " << seeds_k.size() << std::endl;
		file << "hits k size is " << hits_k.size() << std::endl;

//		std::cout << "hits k size is " << hits_k.size() << std::endl;

		if (seeds_k.size() == 0)
			iterate = false;
		if (hits_k.size() < cuts::nseed_hits)
			iterate = false;

		file.close();
	}
	// assign indices
	for (int i=0; i < tracks_k.size(); i++)
        {
                tracks_k[i]->index = i;
        }

}
