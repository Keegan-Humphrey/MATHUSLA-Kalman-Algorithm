#include "globals.hh"
#include "VertexFinder.hh"
#include "statistics.hh"
#include "physics.hh"
#include "LinearAlgebra.hh"
#include "kalman.hh"
#include "kalman-test.hh"
#include <iostream>
#include <fstream>

void VertexFinder::Seed()
{
	for (int n1 = 0; n1 < tracks.size(); n1++)
	{
		for (int n2 = n1 + 1; n2 < tracks.size(); n2++)
		{

			auto tr1 = tracks[n1];
			auto tr2 = tracks[n2];

			if (tr1->closest_approach(tr2) < cuts::seed_closest_approach)
				seeds.push_back(vertex_seed(tr1, tr2));

		} //n2
	}	  //n1

	std::sort(seeds.begin(), seeds.end(), [](vertex_seed a, vertex_seed b) -> bool
			  { return a.score() < b.score(); });
} //VF:Seed

void VertexFinder::Seed_k()
{
	for (int n1 = 0; n1 < tracks_k.size(); n1++)
	{
		for (int n2 = n1 + 1; n2 < tracks_k.size(); n2++)
		{

			auto tr1 = tracks_k[n1];
			auto tr2 = tracks_k[n2];

			if (tr1->closest_approach(tr2) < cuts::seed_closest_approach)
				seeds_k.push_back(vertex_seed(tr1, tr2));

		} //n2
	}	  //n1

	std::sort(seeds_k.begin(), seeds_k.end(), [](vertex_seed a, vertex_seed b) -> bool
			  { return a.score() < b.score(); });
}


void VertexFinder::Seed_k_m()
{
	std::ofstream file;
        file.open("print.txt", std::ios_base::app);

	for (int n1 = 0; n1 < tracks_k_m.size(); n1++)
	{
		for (int n2 = n1 + 1; n2 < tracks_k_m.size(); n2++)
		{

			auto tr1 = tracks_k_m[n1];
			auto tr2 = tracks_k_m[n2];

			if (tr1->closest_approach(tr2) < cuts::seed_closest_approach)
			{
//				vertex_seed possible_seed = vertex_seed(tr1, tr2);
//				double y_seed = possible_seed.guess_k()[1];

//				file << "seed " << y_seed << " tracks: " << tr1->y0 << ", " << tr2->y0 << std::endl;

//				if (tr1->y0 > y_seed && tr2->y0 > y_seed) // check if right topology
//					seeds_k_m.push_back(possible_seed);
					seeds_k_m.push_back(vertex_seed(tr1, tr2));
			}
		} //n2
	}	  //n1

	std::sort(seeds_k_m.begin(), seeds_k_m.end(), [](vertex_seed a, vertex_seed b) -> bool
			  { return a.score() < b.score(); });

	file.close();
} //VF:Seed

void VertexFinder::FindVertices()
{

	if (seeds.size() < 1)
		return; // no seeds

	while (seeds.size() > 0 and tracks.size() > 0)
	{

		std::vector<physics::track *> used_tracks = {};
		std::vector<physics::track *> unused_tracks = {};

		auto current_seed = seeds[0];

		seeds.erase(seeds.begin());

		for (auto tr : tracks)
		{
			if (current_seed.closest_approach(tr) < cuts::closest_approach_add)
			{
				used_tracks.push_back(tr);
			}
			else
			{
				unused_tracks.push_back(tr);
			}
		}

		if (used_tracks.size() < 2)
			continue;

		VertexFitter fitter;
		auto status = fitter.fit(used_tracks, current_seed.guess().std());

		if (status == false or fitter.merit() > cuts::vertex_chi2)
		{
			if (status == false)
				noConverge += 1;
			if (fitter.merit() > cuts::vertex_chi2)
				missedChi2 += 1;
			continue;
		}
		double cos_opening_angle = -1.0;
		if (used_tracks.size() == 2)
		{

			auto tr1 = used_tracks[0];
			auto tr2 = used_tracks[1];

			cos_opening_angle = tr1->vx * tr2->vx + tr1->vy * tr2->vy + tr1->vz * tr2->vz;
			cos_opening_angle = cos_opening_angle / (tr1->beta() * tr2->beta() * constants::c * constants::c);
		}

		auto good_vertex = new physics::vertex(fitter.parameters, cos_opening_angle);


		good_vertex->CovMatrix(fitter.cov_matrix, fitter.npar);
		good_vertex->merit(fitter.merit());
		for (auto tr : used_tracks) good_vertex->track_indices.push_back(tr->index);

		vertices.push_back(good_vertex);
		tracks = unused_tracks;
	}
}

void VertexFinder::FindVertices_k_hybrid()
{

	if (seeds_k.size() < 1)
		return; // no seeds

	while (seeds_k.size() > 0 and tracks_k.size() > 0)
	{

		std::vector<physics::track *> used_tracks = {};
		std::vector<physics::track *> unused_tracks = {};

		auto current_seed = seeds_k[0];

		seeds_k.erase(seeds_k.begin());

		for (auto tr : tracks_k)
		{
			if (current_seed.closest_approach(tr) < cuts::closest_approach_add)
			{
				used_tracks.push_back(tr);
			}
			else
			{
				unused_tracks.push_back(tr);
			}
		}

		if (used_tracks.size() < 2)
			continue;

		VertexFitter fitter;
		auto status = fitter.fit(used_tracks, current_seed.guess().std());

		if (status == false or fitter.merit() > cuts::vertex_chi2)
		{
			if (status == false)
				noConverge += 1;
			if (fitter.merit() > cuts::vertex_chi2)
				missedChi2 += 1;
			continue;
		}
		double cos_opening_angle = -1.0;
		if (used_tracks.size() == 2)
		{

			auto tr1 = used_tracks[0];
			auto tr2 = used_tracks[1];

			cos_opening_angle = tr1->vx * tr2->vx + tr1->vy * tr2->vy + tr1->vz * tr2->vz;
			cos_opening_angle = cos_opening_angle / (tr1->beta() * tr2->beta() * constants::c * constants::c);
		}

		auto good_vertex = new physics::vertex(fitter.parameters, cos_opening_angle);

		for (auto track : used_tracks)
		{
			good_vertex->track_indices.push_back(track->index);
		}

		//good_vertex->CovMatrix(fitter.cov_matrix, fitter.npar);
		good_vertex->merit(fitter.merit());
		vertices_k.push_back(good_vertex);
		tracks_k = unused_tracks;
	}
}

void VertexFinder::FindVertices_k_m_hybrid()
{	if (seeds_k_m.size() < 1)
		return; // no seeds

	while (seeds_k_m.size() > 0 and tracks_k_m.size() > 0)
	{

		std::vector<physics::track *> used_tracks = {};
		std::vector<physics::track *> unused_tracks = {};

		auto current_seed = seeds_k_m[0];

		seeds_k_m.erase(seeds_k_m.begin());

		for (auto tr : tracks_k_m)
		{
			if (current_seed.closest_approach(tr) < cuts::closest_approach_add)
			{
				used_tracks.push_back(tr);
			}
			else
			{
				unused_tracks.push_back(tr);
			}
		}
//		std::cout << "test num tracks" << std::endl;

		if (used_tracks.size() < 2)
			continue;

		VertexFitter fitter;
		auto status = fitter.fit(used_tracks, current_seed.guess().std());

		if (status == false or fitter.merit() > cuts::vertex_chi2)
		{
			if (status == false)
				noConverge += 1;
			if (fitter.merit() > cuts::vertex_chi2)
				missedChi2 += 1;
			continue;
		}

//		std::cout << "test track fit" << std::endl;

		double cos_opening_angle = -1.0;
		if (used_tracks.size() == 2)
		{

			auto tr1 = used_tracks[0];
			auto tr2 = used_tracks[1];

			cos_opening_angle = tr1->vx * tr2->vx + tr1->vy * tr2->vy + tr1->vz * tr2->vz;
			cos_opening_angle = cos_opening_angle / (tr1->beta() * tr2->beta() * constants::c * constants::c);
		}

		auto good_vertex = new physics::vertex(fitter.parameters, cos_opening_angle);

		for (auto track : used_tracks)
		{
			good_vertex->track_indices.push_back(track->index);
		}

		good_vertex->CovMatrix(fitter.cov_matrix, fitter.npar);
		good_vertex->merit(fitter.merit());
		vertices_k_m.push_back(good_vertex);
		tracks_k_m = unused_tracks;

//		std::cout << "vertex_made it " << std::endl;
	}
}
void VertexFinder::FindVertices_k()
{
	std::ofstream file;
	file.open("print.txt", std::ios_base::app);

	file << "-------------------- Begin Vertexing -------------" << std::endl;

	file.close();

	//std::cout << "seeds left m" << seeds_k_m.size() << std::endl;
	//std::cout << "tracks left m" << tracks_k_m.size() << std::endl;
//	std::cout << "seeds left " << seeds_k.size() << std::endl;
//	std::cout << "tracks left " << tracks_k.size() << std::endl;

//	if (seeds_k.size() < 1)
//		return; // no seeds_k
	if (seeds_k_m.size() < 1)
		return;

	for (auto tr : tracks_k_m) file << "track indices: " << tr->index << std::endl;

	while (seeds_k_m.size() > 0 and tracks_k_m.size() > 0)
//	while (seeds_k.size() > 0 and tracks_k.size() > 0)
	{
		//std::cout << "seeds left " << seeds_k.size() << std::endl;
		//std::cout << "tracks left " << tracks_k.size() << std::endl;
//		std::cout << "seeds left " << seeds_k_m.size() << std::endl;

		file.open("print.txt", std::ios_base::app);

		file << "-------------------- New Vertex Seed -------------" << std::endl;

		file.close();

		std::vector<physics::track *> used_tracks = {};
		std::vector<physics::track *> unused_tracks = {};

		auto current_seed = seeds_k_m[0];
//		auto current_seed = seeds_k[0];

		seeds_k_m.erase(seeds_k_m.begin());

//		seeds_k.erase(seeds_k.begin());

//		return;

//		std::cout << "Test 1" << std::endl;

//		unused_tracks = tracks_k_m;
		used_tracks = tracks_k_m;

		double drops = -1;

		int i = 0;

		while (drops != 0)
		{
		file.open("print.txt", std::ios_base::app);

		file << "track indices ";
		for (auto tr : tracks_k_m) file << tr->index << ", ";
		file << std::endl;

		if (i > 10) std::cout << i << std::endl;
		file << "Iteration: " << i << std::endl;

		file << "Number of tracks left " << unused_tracks.size() << std::endl;

		file.close();

		kalman_vertex kfv;
		kfv.dropping = true;
		kfv.vertexer(used_tracks, &current_seed);
//		kfv.vertexer(unused_tracks, &current_seed);
//		kfv.vertexer(tracks_k_m, &current_seed);
//		kfv.vertexer(tracks_k, &current_seed);

//		std::cout << "Test 2" << std::endl;

		//return;
		//std::cout << "done fitter " << std::endl;

		//if (kfv.added_tracks.size() == 0) continue;

		file.open("print.txt", std::ios_base::app);

		if (kfv.status != 2) {
			file << "Vertex failed during fitting! status " << kfv.status << std::endl;
			break;
			//continue;
		}

//		unused_tracks = kfv.unadded_tracks; // just commented out ~~~~~~~~~~~
//		unused_tracks.insert(unused_tracks.end(),kfv_2.unadded_tracks.begin(),kfv_2.unadded_tracks.end());

//		used_tracks = kfv.added_tracks; // just commented out~~~~~~~~~~~~
		used_tracks = {}; // clear used tracks, only push back good ones

		drops = 0;

		// drop tracks that don't make the beta cut
		for (int k=0; k < kfv.added_tracks.size(); k++) {
//			double v = 0;
//			for (int i=0; i < 3; i++) v += std::pow(kfv.added_tracks[k]->q_s[i],2);
//			v = std::sqrt(v);

//			if (std::abs(kfv.pulls_v[k]) < kalman::pull_cut_drop) {

			double beta = (kfv.added_tracks[k]->q_s).norm() / constants::c;

			if (!(kalman::v_cut_drop[0] < beta && beta < kalman::v_cut_drop[1])) {
//			if (!(kalman::v_cut_drop[0] < kfv.pulls_v[k] && kfv.pulls_v[k] < kalman::v_cut_drop[1])) {
				file << "Track dropped with beta " << beta << std::endl;
				unused_tracks.push_back(kfv.added_tracks[k]);
				drops++;
			}
			else {
				used_tracks.push_back(kfv.added_tracks[k]);
			}
		}
/**/

		file << "drops is " << drops << std::endl;

		file << "unused track indices ";
		for (auto tr : unused_tracks) file << tr->index << ", ";
		file << std::endl;

		file << "used track indices ";
		for (auto tr : used_tracks) file << tr->index << ", ";
		file << std::endl;

		if (i > 10) std::cout << "tracks " << used_tracks.size() << std::endl;
		file << "tracks " << used_tracks.size() << std::endl;

		if (used_tracks.size() < 2) {
			file << "Sorry not enough tracks" << std::endl;
			break;
			//continue;
		}

		file.close();

		//if (failed) break;

		if (i == 25) {std::cout << "i break" << std::endl; break;}
		i++;

		if (drops != 0) continue; // only evaluate vertex if no tracks were dropped, otherwise keep going

		//kalman_vertex kfv_2;
		//kfv_2.dropping = false;
		//kfv_2.vertexer(used_tracks, &current_seed);

		file.open("print.txt", std::ios_base::app);

		/*
		if (kfv.status != 2) {
			file << "Vertex failed during fitting!" << std::endl;
			continue;
		}
		*/

//		used_tracks = kfv_2.added_tracks;
//		unused_tracks.insert(unused_tracks.end(),kfv_2.unadded_tracks.begin(),kfv_2.unadded_tracks.end());

/*		if (used_tracks.size() < 2) {
			file << "Sorry not enough tracks" << std::endl;
			break;
//			continue;
		}*/

		// start of used to be kfv_2

		int n_track_params = 4;
        	int ndof = ((4.0 + 3.0) * kfv.added_tracks.size() - n_track_params ); // 4 x_k and 3 q_k parameters
	        if (ndof < 1) ndof = 1;

		if (kfv.chi_v / ndof > cuts::kalman_vertex_chi) {
			file << "Vertex failed with chi " << kfv.chi_v / ndof << "!" << std::endl;
			break;
//			continue;
		}

		file << "Vertex made it with chi " << kfv.chi_v / ndof << "!" << std::endl;
		file.close();


		auto good_vertex = new physics::vertex(kfv.x_s, -2.0);

                for (auto track : used_tracks)
                {
                        good_vertex->track_indices.push_back(track->index);

			good_vertex->q_f.push_back(track->q_f);
			good_vertex->D_f.push_back(track->D_f);

			good_vertex->q_s.push_back(track->q_s);
			good_vertex->D_s.push_back(track->D_s);
                }

                //good_vertex->CovMatrix(fitter.cov_matrix, fitter.npar);
                //good_vertex->merit(fitter.merit());

		vertices_k_m.push_back(good_vertex);

                tracks_k_m = unused_tracks;

		} // drop while loop

//		std::cout << "Test 3" << std::endl;

		//return;
	}
	//std::cout << " done vertexing" << std::endl;

}

std::vector<physics::track *> VertexFitter::track_list = {};
std::vector<double> VertexFitter::parameters = {};
std::vector<double> VertexFitter::parameter_errors = {};
double VertexFitter::_merit = 0.;
double VertexFitter::cov_matrix[VertexFitter::npar][VertexFitter::npar];

bool VertexFitter::bad_fit = false;
void VertexFitter::nll(int &npar, double *gin, double &f, double *pars, int iflag)
{
	std::ofstream file;
	file.open("print.txt", std::ios_base::app);

	using Vector = vector::Vector;
	double _x = pars[0];
	double _y = pars[1];
	double _z = pars[2];
	double _t = pars[3];

	double error = 0.0;

	int n = 0;

	for (auto track : VertexFitter::track_list)
	{

		double dist = track->distance_to(Vector(_x, _y, _z), _t);

		double err = track->err_distance_to(Vector(_x, _y, _z), _t);

		error += 0.5 * (dist / err) * (dist / err) + TMath::Log(err);

		if (isnan(error))
		{
			bad_fit = true;
			file << " Bad fit! " << std::endl;
			file << dist << " " << err << std::endl;
//			track->CovMatrix().Print();
			return;
		}
	}

	f = error;
	file.close();
}
