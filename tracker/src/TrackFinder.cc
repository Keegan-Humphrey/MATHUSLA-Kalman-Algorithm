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
#include "par_handler.hh"
#include "Math/ProbFunc.h"

void TrackFinder::Seed()
{
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

			double ds_2 = c_score(hits[first], hits[second]);

			//double dr = (hits[first].PosVector() - hits[second].PosVector()).Magnitude() / constants::c; // [ns]

			if (ds_2 > par_handler->par_map["seed_interval"])
				continue;

			seeds.push_back(seed(hits[first], hits[second]));
			seeds_k.push_back(seed(hits[first], hits[second]));

		} //"second" loop
	}	  //"first" loop

	score_seeds();

} //TF::Seed

void TrackFinder::FindTracks()
{
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

		int min_index = min_seed();
		auto current_seed = seeds[min_index];

		seeds.erase(seeds.begin() + min_index); //delete the seed so that it isn't used again

		std::vector<physics::digi_hit *> track_pts;
		std::vector<physics::digi_hit *> unused_hits;

		//now we use the seed::residual function to get the residual from the seed of all of the available hits.
		//if the residual is good, we add it to the track_pts vector

		for (auto hit : hits)
		{
			if ((current_seed.timeless_residual(hit) < cuts::seed_residual or current_seed.distance_to_hit(hit) < cuts::distance_to_hit) and current_seed.time_difference(hit) < cuts::seed_time_difference)
			{
				track_pts.push_back(hit);
			}
			else
			{

				unused_hits.push_back(hit);
			}
		}

		//we check that there are at least cuts::nseed_hits hits in the track_pts
		//if not, we move on to the next seed

		if (track_pts.size() < cuts::nseed_hits)
		{
			continue;
		}

		//at this point, all we have done is erase the seed, and none of the hits have been modified
		TrackFitter fitter;
		auto track_status = fitter.fit(track_pts, current_seed.guess());

		if (track_status == false) {
			file.close();
			continue; //fit failed
		}

		auto current_track = new physics::track(fitter.parameters, fitter.parameter_errors);
		for (auto hit : track_pts)
			current_track->AddHit(hit);

		//first pass track is now cconstructed

		std::vector<physics::digi_hit *> second_unused_hits;

		for (auto hit : unused_hits)
		{
			if ((current_track->untimed_residual(hit) < cuts::residual_add or current_track->distance_to_hit(hit) < cuts::distance_to_hit) and current_track->time_difference(hit) < cuts::seed_time_difference)
			{
				current_track->AddHit(hit);
			}
			else
			{
				second_unused_hits.push_back(hit);
			}
		}

		//NOW WE REFIT AGAIN

		fitter.fit(current_track->hits, current_track->parameters());
		current_track->parameters(fitter.parameters);
		current_track->par_errors(fitter.parameter_errors);

		std::vector<physics::digi_hit *> good_hits;
		for (auto hit : current_track->hits)
		{

			if (current_track->untimed_residual(hit) > cuts::residual_drop or current_track->time_difference(hit) > cuts::time_difference_drop)
			{
				second_unused_hits.push_back(hit);
			}
			else
			{
				good_hits.push_back(hit);
			}
		}

		if (good_hits.size() < cuts::nseed_hits)
		{
			continue;
		}

		current_track->hits = good_hits;

		fitter.fit(current_track->hits, current_track->parameters());
		current_track->parameters(fitter.parameters);
		current_track->par_errors(fitter.parameter_errors);
		current_track->CovMatrix(fitter.cov_matrix, fitter.npar);

		if (current_track->nlayers() >= cuts::track_nlayers and current_track->chi2_per_dof() < cuts::track_chi2 and current_track->hits.size() > cuts::ntrack_hits)
		{
			tracks.push_back(current_track);
		}
		else
		{
			delete current_track;
			continue;
		}

		hits = second_unused_hits;

		if (seeds.size() == 0)
			iterate = false;
		if (hits.size() < cuts::nseed_hits)
			iterate = false;

	} //while iterate

	//MergeTracks();
	//CleanTracks();

	int total_tracked_points = 0;

	for (auto track : tracks)
		total_tracked_points += track->hits.size();

	if ((total_hits - hits.size() - total_tracked_points) != 0)
	{
		std::cout << "total: " << total_hits << std::endl;
		std::cout << "used: " << total_tracked_points << std::endl;
		std::cout << "unused:" << hits.size() << std::endl;
	}

	if (j++ > MAX_ITS)
		iterate = false;

}

void TrackFinder::MergeTracks()
{
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

//			if (distance > cuts::merge_distance or cos_theta < cuts::merge_cos_theta)
			if (distance > par_handler->par_map["merge_distance"] or cos_theta < par_handler->par_map["merge_cos_theta"])
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
}

void TrackFinder::MergeTracks_k()
{
	//at this point, all of the points have been fit to tracks. At this point, we will perform the track merging step.

	if (tracks_k_m.size() == 0 or tracks_k_m.size() == 1)
		return;

	std::vector<int> deleted_tracks = {};

	Geometry geo;

	for (int first_track = 0; first_track < tracks_k_m.size(); first_track++)
	{

		for (int second_track = first_track + 1; second_track < tracks_k_m.size(); second_track++)
		{
			auto tr1 = tracks_k_m[first_track];
			auto tr2 = tracks_k_m[second_track];

			bool mergebool = false;

			auto first_layer = tr1->layers()[0] < tr2->layers()[0] ? tr1->layers()[0] : tr2->layers()[0];

			auto d1 = tr1->direction();
			auto d2 = tr2->direction();

			int tr1_missing_hits = 0;
			int tr2_missing_hits = 0;

			auto cos_theta = d1 ^ d2;

			double distance = tr1->closest_approach(tr2);

//			if (distance < cuts::merge_distance and cos_theta > cuts::merge_cos_theta){
			if (distance < par_handler->par_map["merge_distance"] and cos_theta > par_handler->par_map["merge_cos_theta"]){
				mergebool = true;
			}

			if ((tr1->_holes.size() >= 3 or tr2->_holes.size() >= 3) and (distance < 150 and cos_theta > 0.95)){
				//std::cout << "first holes" << std::endl;
				mergebool = true;
			}

			if ((tr1->_holes.size() >= 2 or tr2->_holes.size() >= 2) and (distance < 100 and cos_theta > 0.99)){
				//std::cout << "second holes" << std::endl;
				mergebool = true;
			}
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

			//if (!merge) // only affects above condition
			//	continue;
			if (!mergebool){continue;} // from above loop

			//std::cout << "made it passed mergebool" << std::endl;

			Vector CA_mid = tr1->closest_approach_midpoint(tr2);
			if (!geo.inBox(CA_mid.x, CA_mid.y, CA_mid.z)) {
				std::cout << "CA is outside of detector" << std::endl;
				continue;
			}

			//std::cout << "made it passed CA" << std::endl;

			seed artificial_seed = seed(tr1->hits[0], tr1->hits[1]);

			// we've decided that the tracks should be merged
			// replace the tracks by a track made of a fit of their combined hits
			//
			// *** now we only take the hits from the first track
			// since the tracks are similar we can just use one track
			// *** the refit can eventually be taken out all together
			//for (auto hit : tr2->hits)
			//	tr1->AddHit(hit);

			kalman_track kft;
			kft.par_handler = par_handler;
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

	tracks_k_m = good_tracks;

	//at this point, the list of tracks is finalized. Now they will be indexed:
	for (int trackn=0; trackn<tracks_k_m.size(); trackn++)
	{
		tracks_k_m[trackn]->index = trackn;
	}
	file.close();
}

void TrackFinder::CleanTracks()
{

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

		} //second track
	}	  //first track
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
					missing_layers.push_back(expected_index);
				}
			} //if missing
		}
		track->set_holes(missing_layers);
	}
}

void TrackFinder::CalculateMissingHits(Geometry *geo)
{

	for (auto track : tracks)
	{
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
					missing_layers.push_back(expected_index);
				}
			} //if missing
		}
		track->missing_layers(missing_layers);
	}
	//calculate for kalman tracks

	for (auto track : tracks_k_m)
	{
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
					missing_layers.push_back(expected_index);
				}
			} //if missing
		}

		track->missing_layers(missing_layers);
	}
}

void TrackFinder::FindTracks_kalman()
{

	if (par_handler->par_map["debug"] == 1) std::cout << "Number of seeds:" <<  seeds_k.size() << std::endl;

	if (seeds_k.size() == 0)
		return; //no seeds found in initial seeding, will be retried with <c travel

	int index = 0;
	int total_hits = hits_k.size();
	bool iterate = true;
	int j = 0;
	int MAX_ITS = 25;

	Stat_Funcs sts;

	while (iterate)
	{
		if (seeds_k.size() == 0)
			break;
		if (hits_k.size() == 0)
			break;

		int min_index = min_seed_k();
		auto current_seed = seeds_k[min_index];

		seeds_k.erase(seeds_k.begin() + min_index); //delete the seed so that it isn't used again

		if (par_handler->par_map["debug"] == 1) std::cout << "New Seed ---------" << std::endl;

		// check if the first seed hit is in the hit pool
		bool used = !seed_unused(current_seed); // double negative! used = not unused

		clear_vecs();

		std::vector<int> king_move_inds;

		// Construct the first filter (to find good hits from seed)
		kalman_track kf_find;
		kf_find.par_handler = par_handler;
		kf_find.finding = true;
		kf_find.dropping = true;
		kf_find.seed_was_used = used;

		if (par_handler->par_map["debug"] == 1) std::cout << "first fit" << std::endl;

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

		undropped_hits = kf_find.found_hits;
		unused_hits = kf_find.unadded_hits;

//		int drops = -1;
		int i = 0;

		bool failed = false;

//		while (drops != 0)
//		{
			kalman_track kft_;
			kft_.par_handler = par_handler;
        	        kft_.finding = false;
                	kft_.dropping = true;
	                kft_.seed_was_used = used;
			kft_.unadded_hits = unused_hits;

			if (par_handler->par_map["debug"] == 1) std::cout << "second fit" << std::endl;

                	kft_.kalman_all(undropped_hits, &current_seed);

			if (kft_.status == -1)
        	        {
				failed = true;
	                        continue;
        	        }

			good_hits = kft_.added_hits;

			undropped_hits.clear();

			if (kft_.status == -2)
			{
				failed = true;
				continue;
			}

			if (good_hits.size() < cuts::track_nlayers)
			{
				failure_reason[4] += 1;
			}

			unused_hits = kft_.unadded_hits;
			unused_hits.insert(unused_hits.end(), kf_find.unadded_hits.begin(), kf_find.unadded_hits.end());

//			int drops = 0;

			double ndof = good_hits.size();
			ndof = ndof > 1.0 ? 4.0 * ndof - 6.0 : 1.0;

			//dropping hits
			for (int n = 0; n < good_hits.size(); n++)
			{
				// make an eigenvector for the velocity at the lowest (in y) state of the track
				Eigen::VectorXd v(3);
				v << kft_.v_s_list[n][0], kft_.v_s_list[n][1], kft_.v_s_list[n][2];

//				if (kft_.chi_s[n] > cuts::kalman_chi_s
				//if (sts.chi_prob_eld(kft_.chi_s[n],ndof) > par_handler->par_map["kalman_pval_s"]
//				if (kft_.chi_s[n] > par_handler->par_map["kalman_chi_s"]
				if (ROOT::Math::chisquared_cdf(kft_.chi_s[n], ndof) >= par_handler->par_map["kalman_pval_drop"]
			           || !(par_handler->par_map["kalman_v_drop[0]"] < v.norm() / constants::c
				   && v.norm() / constants::c < par_handler->par_map["kalman_v_drop[1]"]))
//			           || !(cuts::kalman_v_drop[0] < v.norm() / constants::c && v.norm() / constants::c < cuts::kalman_v_drop[1]))
				{
					if (par_handler->par_map["debug"] == 1) {
						std::cout << "hit dropped with y " << good_hits[n]->y << " chi " << ROOT::Math::chisquared_cdf(kft_.chi_s[n], ndof) <<
							" v " << v.norm() / constants::c << std::endl;
					}
					unused_hits.push_back(good_hits[n]);
//					drops++;
				}
				else
				{
					undropped_hits.push_back(good_hits[n]);
				}
			}

			if (undropped_hits.size() < cuts::track_nlayers)
			{
				failed = true;
//				break;
				continue;
			}

		//} // dropping while loop

		// Do a final hit with bad hits removed
		kalman_track kft_2;
		kft_2.par_handler = par_handler;
		kft_2.finding = false;
		kft_2.dropping = false;

		if (par_handler->par_map["debug"] == 1) std::cout << "third fit" << std::endl;

		kft_2.kalman_all(undropped_hits, &current_seed);

		if (kft_2.status != 2)
		{
			continue;
		}

		std::vector<double> diag_cov;
		for (int i = 0; i < 7; i++) diag_cov.push_back(kft_2.track_cov[i][i]);

		//auto current_track = new physics::track(kft_2.x_s, kft_2.P_s);
		auto current_track = new physics::track(kft_2.x_s, diag_cov);
		current_track->hits = undropped_hits;

		current_track->x_scats = kft_2.x_scat;
		current_track->z_scats = kft_2.z_scat;

		current_track->chi_f = kft_2.chi_f;
		current_track->chi_s = kft_2.chi_s;

		current_track->estimate_list = kft_2.x_s_list;
		current_track->P_s = kft_2.P_s0;

		current_track->king_move_inds = kft_2.king_move_inds;

		/*
		static double cov_matrix[7][7];
		for (int i = 0; i < 7; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				//cov_matrix[i][j] = kft_2.track_cov[i][i];
				//if (i == j)
				//{
				//	cov_matrix[i][j] = std::pow(kft_2.P_s[i],2);
				//}
				else
				{
					cov_matrix[i][j] = 0;
				}

			}
		}
		*/
		//current_track->CovMatrix(cov_matrix, 7);

		//static double cov_matrix[7][7];
		//cov_matrix = kft_2.track_cov;
		current_track->CovMatrix(kft_2.track_cov, 7);

		//double ndof = kft_2.chi_s.size();
		ndof = kft_2.chi_s.size();
		ndof = ndof > 1.0 ? 4.0 * ndof - 6.0 : 1.0;

		for (auto chi : kft_2.chi_f) {
			local_chi_f.push_back(chi);

			//double p_val = sts.chi_prob(chi,ndof);
                        //local_chi_f.push_back(p_val);
		}
		local_chi_f.push_back(-1);

		for (auto chi : kft_2.chi_s) {
			local_chi_s.push_back(chi);

			//double p_val = sts.chi_prob(chi,ndof);
			//local_chi_s.push_back(p_val);
		}
		local_chi_s.push_back(-1);

		// calculate chi per ndof from sum of chi increments in the track
		double chi_sum = 0;
		for (auto chi : kft_2.chi_s)
			chi_sum += chi;

		//chi_sum = chi_sum / (4.0 * kft_2.chi_s.size() - 6.0); // always non_negative due to num_hit cut (if track can be passed)
		chi_sum = chi_sum / ndof;

//		if (current_track->nlayers() >= cuts::track_nlayers && chi_sum < cuts::kalman_track_chi)
		if (current_track->nlayers() >= cuts::track_nlayers && chi_sum < par_handler->par_map["kalman_track_chi"])
		{
			tracks_k.push_back(current_track);

			if (par_handler->par_map["debug"] == 1) std::cout << "Track made it" << std::endl;
		}
		else
		{
			delete current_track;
			continue;
		}


		if (failed) continue;

		hits_k = unused_hits;

		if (seeds_k.size() == 0)
			iterate = false;
		if (hits_k.size() < cuts::nseed_hits)
			iterate = false;

	}
	// assign indices
	for (int i=0; i < tracks_k.size(); i++)
        {
                tracks_k[i]->index = i;
        }

}
