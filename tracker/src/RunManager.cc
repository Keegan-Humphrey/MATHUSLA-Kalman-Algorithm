#include <iostream>
#include "RunManager.hh"
#include "TreeHandler.hh"
#include "TrackFinder.hh"
#include "Digitizer.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "NoiseMaker.hh"
#include "par_handler.hh"

int RunManager::StartTracking()
{
	TreeHandler _handler(_InputTree_Name, _InputFile_Name, _OutputTree_Name, OutFileName());
	if (_handler.IsNull()) {
		std::cout << "Sorry, I couldn't open that file" << std::endl;
		return 0;
	} 

	TH = &_handler;
	int events_handled = 0;

	int made_its = 0;
	int made_its_k = 0;

	int events_w_tracks = 0;
	int events_w_k_tracks = 0;

	int neg_covs = 0;

	int verts = 0;
	int verts_k = 0;
	int verts_k_m = 0;

	int dropped_hits = 0;
	int floor_wall_hits = 0;

	std::vector<int> zeros(9, 0);
	std::vector<int> failures_k = zeros;

	ParHandler hndlr;
	hndlr.Handle();

	if (hndlr.par_map["branch"] == 1.0) std::cout << "Running in Cosmic Mode" << std::endl;

	_digitizer->par_handler = &hndlr;
	_tracker->par_handler = &hndlr;
	_vertexer->par_handler = &hndlr;

    NoiseMaker::preDigitizer();
    
	_digitizer->InitGenerators();

	while (TH->Next() >= 0)
	{
		if (events_handled >= hndlr.par_map["end_ev"]) //cuts::end_ev)
		{
			break;
		}
		if (events_handled >= hndlr.par_map["start_ev"]) //cuts::start_ev)
		{

			if ((events_handled - 1) % 1000 == 0 || hndlr.par_map["debug"] == 1)
				std::cout << "Event is " << events_handled - 1 << std::endl;

			TotalEventsProcessed++;
			_digitizer->clear();
			_tracker->clear();
			_vertexer->clear();

			// copying the data to the new tree, and loading all the variables, incrementing index
			TH->LoadEvent();

			//adding all hits of the tree into the digitizer
			for (int n_hit = 0; n_hit < TH->sim_numhits; n_hit++)
			{
				physics::sim_hit *current = new physics::sim_hit(TH, n_hit);
				if (hndlr.par_map["branch"] == 1.0) {
					current->x += detector::COSMIC_SHIFT[0];
					current->y += detector::COSMIC_SHIFT[1];
					current->z += detector::COSMIC_SHIFT[2];
				}
				_digitizer->AddHit(current);
			}

			_digitizer->ev_num = events_handled; // used to vary the seed, otherwise for events < 1 s seperatred in run time will have the same seed
                                            // ^^ not true anymore, can likely be removed?
			std::vector<physics::digi_hit *> digi_list = _digitizer->Digitize();

			TH->ExportDigis(digi_list, _digitizer->seed);

			// digis now finished and stored in tree!!!
			// now, we begin the seeding algorithm

			// remove this carefully in TrackFinder.cc
			_tracker->failure_reason = zeros;

			_tracker->hits = digi_list;
			_tracker->hits_k = digi_list;
			_tracker->Seed();
			//_tracker->FindTracks();
			_tracker->FindTracks_kalman();

			// copy kalman tracks for merging
			for (auto t : _tracker->tracks_k)
			{
				physics::track *temp = new physics::track(*t);
				_tracker->tracks_k_m.push_back(temp);
			}

			_tracker->CalculateMissingHits(_digitizer->_geometry);
			_tracker->MergeTracks();

			if (hndlr.par_map["merge_cos_theta"] != -2) {
				_tracker->CalculateHoles(_digitizer->_geometry);
				_tracker->MergeTracks_k();
			}

			_tracker->CalculateMissingHits(_digitizer->_geometry);

			made_its += _tracker->tracks.size();
			made_its_k += _tracker->tracks_k_m.size();

			TH->ExportTracks(_tracker->tracks);
			//TH->ExportTracks_k(_tracker->tracks_k);
			TH->ExportTracks_k_m(_tracker->tracks_k_m);
			//TH->Export_kalman_info(_tracker->local_chi_f, _tracker->local_chi_s);

			_vertexer->tracks = _tracker->tracks;
			_vertexer->tracks_k = _tracker->tracks_k;
			_vertexer->tracks_k_m = _tracker->tracks_k_m;

			_vertexer->Seed();
			//_vertexer->Seed_k();
			_vertexer->Seed_k_m();
			_vertexer->FindVertices();
			_vertexer->FindVertices_k_m_hybrid();

			//_vertexer->FindVertices_k();

			verts += _vertexer->vertices.size();
			//verts_k += _vertexer->vertices_k.size();
			verts_k_m += _vertexer->vertices_k_m.size();

			TH->ExportVertices(_vertexer->vertices);
			//TH->ExportVertices_k(_vertexer->vertices_k);
			TH->ExportVertices_k_m(_vertexer->vertices_k_m);

			TH->Fill();

			if (_tracker->tracks.size() > 0) events_w_tracks++;
			if (_tracker->tracks_k_m.size() > 0) events_w_k_tracks++;

			dropped_hits += _digitizer->dropped_hits;
			floor_wall_hits += _digitizer->floor_wall_hits;

		}

		events_handled++;
	}

	if (hndlr.file_opened) {
		std::cout << made_its << " Linear tracks made it" << std::endl;
		std::cout << made_its_k << " Kalman tracks made it" << std::endl;
		std::cout << verts << " Linear vertices made it" << std::endl;
		//std::cout << verts_k << " Kalman vertices made it" << std::endl;
		std::cout << verts_k_m << " Merged Kalman vertices made it" << std::endl;
		std::cout << events_w_tracks << " Events had a linear track" << std::endl;
		std::cout << events_w_k_tracks << " Events had a kalman track" << std::endl;

		TH->Write();

		std::cout << "Tracked " << TotalEventsProcessed << " Events" << std::endl;

//		std::cout << "# Dropped Hits / # Floor Wall Hits = " << (double) dropped_hits / floor_wall_hits;

	}

	return 0;
}
