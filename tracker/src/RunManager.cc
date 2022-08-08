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

	int made_its_k = 0;
	int made_its_c_b = 0;

	int events_w_k_tracks = 0;
	int events_c_b_tracks = 0;

	int verts_k_m = 0;
	int verts_c_b = 0;
	//TODO: should I copy these 4 below?
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

			_digitizer->ev_num = events_handled;
			std::vector<physics::digi_hit *> digi_list = _digitizer->Digitize();

			TH->ExportDigis(digi_list, _digitizer->seed);

			// digis now finished and stored in tree!!!
			// now, we begin the seeding algorithm

			// remove this carefully in TrackFinder.cc
			_tracker->failure_reason = zeros;

			_tracker->hits_k = digi_list;
			_tracker->hits_c_b = digi_list; //TODO look into tracker names
			_tracker->Seed();
			_tracker->FindTracks_kalman();
			_tracker->FindTracks_c_b(); //TODO look into tracker names

			// copy kalman tracks for merging
			for (auto t : _tracker->tracks_k)
			{
				physics::track *temp = new physics::track(*t);
				_tracker->tracks_k_m.push_back(temp);
			}

			for (auto t : _tracker->tracks_c_b)
			{
				physics::track *temp = new physics::track(*t);
				_tracker->tracks_c_b.push_back(temp); //TODO look into tracker names
			}
			
			
			if (hndlr.par_map["merge_cos_theta"] != -2) {
				_tracker->CalculateMissingHits(_digitizer->_geometry);
				_tracker->MergeTracks_k();
			}

			if (hndlr.par_map["merge_cos_theta"] != -2) {
				_tracker->CalculateMissingHits_c_b(_digitizer->_geometry);//Tracker name
				_tracker->MergeTracks_c_b();//TODO: Tracker name
			}
			
			_tracker->CalculateMissingHits(_digitizer->_geometry);
			_tracker->CalculateMissingHits_c_b(_digitizer->_geometry); //TODO: tracker name

			made_its_k += _tracker->tracks_k_m.size();
			made_its_c_b += _tracker->tracks_c_b.size(); //TODO: tracker name

			TH->ExportTracks_k_m(_tracker->tracks_k_m);
			TH->ExportTracks_c_b(_tracker->tracks_c_b); //TODO: tracker name

			_vertexer->tracks_k_m = _tracker->tracks_k_m;
			_vertexer->tracks_c_b = _tracker->tracks_c_b; //TODO: vertexer name

			_vertexer->Seed_k_m();
			_vertexer->Seed_c_b(); //TODO: vertexer name
			_vertexer->FindVertices_k_m_hybrid();
			_vertexer->FindVertices_c_b_hybrid(); //TODO: vertexer name

			//_vertexer->FindVertices_k();

			verts_k_m += _vertexer->vertices_k_m.size();
			verts_c_b += _vertexer->vertices_c_b.size(); //TODO: vertexer name

			TH->ExportVertices_k_m(_vertexer->vertices_k_m);
			TH->ExportVertices_c_b(_vertexer->vertices_c_b); //TODO: vertexer name

			TH->Fill();

			if (_tracker->tracks_k_m.size() > 0) events_w_k_tracks++;
			if (_tracker->tracks_c_b.size() > 0) events_c_b_tracks++; //TODO: tracker name

			dropped_hits += _digitizer->dropped_hits;
			floor_wall_hits += _digitizer->floor_wall_hits;

		}

		events_handled++;
	}

	if (hndlr.file_opened) {
		std::cout << made_its_k << " Kalman tracks made it" << std::endl;
		std::cout << verts_k_m << " Merged Kalman vertices made it" << std::endl;
		std::cout << events_w_k_tracks << " Events had a kalman track" << std::endl;
		
		std::cout << made_its_c_b << " Limited Kalman tracks made it" << std::endl;
		std::cout << verts_c_b << " Merged Limited Kalman vertices made it" << std::endl;
		std::cout << events_c_b_tracks << " Events had a limited kalman track" << std::endl;

		TH->Write();

		std::cout << "Tracked " << TotalEventsProcessed << " Events" << std::endl;

//		std::cout << "# Dropped Hits / # Floor Wall Hits = " << (double) dropped_hits / floor_wall_hits;

	}

	return 0;
}
