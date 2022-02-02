#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <Eigen/Dense>

#ifndef TH_DEFINE
#define TH_DEFINE

class TreeHandler{
public:
	//PUT ALL INPUT AND OUTPUT BRANCHES HERE
	TTree* OutputTree;
	TTree* InputTree;
	TFile* OutputFile;
  TFile* InputFile;
	int index = -1;
	int NumEntries;
  bool _Null = false;
  bool IsNull(){return _Null;}

	int Next(){
    index++;
		if (index >= NumEntries) return -1;
		return index;
	}

	int LoadEvent(){
		if (InputTree == nullptr) return -1;
		InputTree->GetEvent(index);
		return 0;
	}

	void Fill(){
		OutputFile->cd();
		gROOT->cd();
		InputTree->GetEvent(index);
		OutputTree->Fill();
	}

	void Write(){
    InputFile->Close();
		OutputFile->cd();
		OutputFile->Write();
		OutputFile->Close();
	}

  void Export_kalman_info(std::vector<double>, std::vector<double>);

  template<class digi_hit>
  void ExportDigis(std::vector<digi_hit*>);

  template<class track>
  void ExportTracks_k(std::vector<track*>);

	template<class track>
  void ExportTracks_k_m(std::vector<track*>);

  template<class track>
  void ExportTracks(std::vector<track*>);

  template<typename vertex>
  void ExportVertices(std::vector<vertex*>);

  template<typename vertex>
  void ExportVertices_k(std::vector<vertex*>);

	template<typename vertex>
  void ExportVertices_k_m(std::vector<vertex*>);

	TreeHandler(TString input_tree_name, TString input_file_name, TString output_tree_name, TString outfile_name)
	{

		InputFile = TFile::Open(input_file_name);
    if (! InputFile) {
      _Null = true;
      return;
    }

		InputTree = (TTree*) InputFile->Get(input_tree_name);

    if (! InputTree){
      _Null = true;
      return;
    }


		InputTree->SetBranchAddress("NumHits", &sim_numhits);
 		InputTree->SetBranchAddress("Hit_energy", &sim_hit_e);
 		InputTree->SetBranchAddress("Hit_time", &sim_hit_t);
 		InputTree->SetBranchAddress("Hit_detId", &sim_hit_detId);
 		InputTree->SetBranchAddress("Hit_particlePdgId", &sim_hit_particlePdgId);
 		InputTree->SetBranchAddress("Hit_G4TrackId", &sim_hit_G4TrackId);
 		InputTree->SetBranchAddress("Hit_G4ParentTrackId", &sim_hit_G4ParentTrackId);
 		InputTree->SetBranchAddress("Hit_x", &sim_hit_x);
 		InputTree->SetBranchAddress("Hit_y", &sim_hit_y);
 		InputTree->SetBranchAddress("Hit_z", &sim_hit_z);
 		InputTree->SetBranchAddress("Hit_particleEnergy", &sim_hit_particleEnergy);
 		InputTree->SetBranchAddress("Hit_particlePx", &sim_hit_px);
 		InputTree->SetBranchAddress("Hit_particlePy", &sim_hit_py);
		InputTree->SetBranchAddress("Hit_particlePz", &sim_hit_pz);
		InputTree->SetBranchAddress("Hit_weight", &sim_hit_weight);
	InputTree->SetBranchAddress("GenParticle_energy", &sim_GenParticle_energy);
	InputTree->SetBranchAddress("GenParticle_pdgid", &sim_GenParticle_pdgid);
	InputTree->SetBranchAddress("GenParticle_index", &sim_GenParticle_index);
	InputTree->SetBranchAddress("GenParticle_G4index", &sim_GenParticle_G4index);
 	InputTree->SetBranchAddress("GenParticle_time", &sim_GenParticle_time);
 	InputTree->SetBranchAddress("GenParticle_x", &sim_GenParticle_x);
 	InputTree->SetBranchAddress("GenParticle_y", &sim_GenParticle_y);
 	InputTree->SetBranchAddress("GenParticle_z", &sim_GenParticle_z);
 	InputTree->SetBranchAddress("GenParticle_px", &sim_GenParticle_px);
 	InputTree->SetBranchAddress("GenParticle_py", &sim_GenParticle_py);
 	InputTree->SetBranchAddress("GenParticle_pz", &sim_GenParticle_pz);
 	InputTree->SetBranchAddress("GenParticle_mass", &sim_GenParticle_mass);


    InputTree->SetBranchStatus("NumGenParticles", 0);
//    InputTree->SetBranchStatus("GenParticle_index", 0);
//    InputTree->SetBranchStatus("GenParticle_G4index", 0);
//    InputTree->SetBranchStatus("GenParticle_pdgid", 0);
    InputTree->SetBranchStatus("GenParticle_status", 0);
//    InputTree->SetBranchStatus("GenParticle_time", 0);
//    InputTree->SetBranchStatus("GenParticle_x", 0);
//    InputTree->SetBranchStatus("GenParticle_y", 0);
//    InputTree->SetBranchStatus("GenParticle_z", 0);
   // InputTree->SetBranchStatus("GenParticle_energy", 0);
//    InputTree->SetBranchStatus("GenParticle_px", 0);
//    InputTree->SetBranchStatus("GenParticle_py", 0);
//    InputTree->SetBranchStatus("GenParticle_pz", 0);

//    InputTree->SetBranchStatus("GenParticle_mass", 0);
    InputTree->SetBranchStatus("GenParticle_pt", 0);
    InputTree->SetBranchStatus("GenParticle_eta", 0);
    InputTree->SetBranchStatus("GenParticle_phi", 0);

    InputTree->SetBranchStatus("GenParticle_mo1", 0);
    InputTree->SetBranchStatus("GenParticle_mo2", 0);
    InputTree->SetBranchStatus("GenParticle_dau1", 0);
    InputTree->SetBranchStatus("GenParticle_dau2", 0);
    InputTree->SetBranchStatus("COSMIC_EVENT_ID", 0);
 		InputTree->SetBranchStatus("COSMIC_CORE_X", 0);
 		InputTree->SetBranchStatus("COSMIC_CORE_Y", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_PRIMARY_ENERGY", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_THETA", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_PHI", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_FIRST_HEIGHT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_ELECTRON_COUNT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_MUON_COUNT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_HADRON_COUNT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_PRIMARY_ID", 0);
 		InputTree->SetBranchAddress("EXTRA_11", &sim_EXTRA_11);
 		InputTree->SetBranchAddress("EXTRA_12", &sim_EXTRA_12);
 		InputTree->SetBranchAddress("EXTRA_13", &sim_EXTRA_13);
 		InputTree->SetBranchAddress("EXTRA_14", &sim_EXTRA_14);
 		InputTree->SetBranchAddress("EXTRA_15", &sim_EXTRA_15);

 		NumEntries = InputTree->GetEntries();

 		OutputFile = new TFile(outfile_name, "RECREATE");
		OutputTree = new TTree(output_tree_name, "MATHUSLA Tree");

		OutputTree->Branch("NumHits", &sim_numhits);
		OutputTree->Branch("Hit_energy", "std::vector<double>", sim_hit_e);
 		OutputTree->Branch("Hit_time", "std::vector<double>", sim_hit_t);
 		OutputTree->Branch("Hit_detId", "std::vector<double>", sim_hit_detId);
 		OutputTree->Branch("Hit_particlePdgId", "std::vector<double>", sim_hit_particlePdgId);
 		OutputTree->Branch("Hit_G4TrackId", "std::vector<double>", sim_hit_G4TrackId);
 		OutputTree->Branch("Hit_G4ParentTrackId", "std::vector<double>", sim_hit_G4ParentTrackId);
 		OutputTree->Branch("Hit_x", "std::vector<double>", sim_hit_x);
 		OutputTree->Branch("Hit_y", "std::vector<double>", sim_hit_y);
 		OutputTree->Branch("Hit_z", "std::vector<double>", sim_hit_z);
 		OutputTree->Branch("Hit_particleEnergy", "std::vector<double>", sim_hit_particleEnergy);
 		OutputTree->Branch("Hit_particlePx", "std::vector<double>", sim_hit_px);
 		OutputTree->Branch("Hit_particlePy", "std::vector<double>", sim_hit_py);
		OutputTree->Branch("Hit_particlePz", "std::vector<double>", sim_hit_pz);
//		OutputTree->Branch("Hit_weight", "std::vector<double>", sim_hit_weight);

      OutputTree->Branch("Digi_numHits", &Digi_numHits);
        OutputTree->Branch("Digi_time", &digi_hit_t);
        OutputTree->Branch("Digi_x", &digi_hit_x);
        OutputTree->Branch("Digi_y", &digi_hit_y);
        OutputTree->Branch("Digi_z", &digi_hit_z);
        OutputTree->Branch("Digi_energy", &digi_hit_e);
        OutputTree->Branch("Digi_px", &digi_hit_px);
        OutputTree->Branch("Digi_py", &digi_hit_py);
        OutputTree->Branch("Digi_pz", &digi_hit_pz);
        OutputTree->Branch("Digi_particle_energy", &digi_particle_energy);
        OutputTree->Branch("Digi_pdg_id", &digi_pdg);
//        OutputTree->Branch("Digi_hitIndices", &digi_hit_indices);

 	//	OutputTree->Branch("NumGenParticles", &sim_NumGenParticles);
 		OutputTree->Branch("GenParticle_index", "std::vector<double>", sim_GenParticle_index);
 		OutputTree->Branch("GenParticle_G4index", "std::vector<double>", sim_GenParticle_G4index);
 		OutputTree->Branch("GenParticle_pdgid", "std::vector<double>", sim_GenParticle_pdgid);
 	//	OutputTree->Branch("GenParticle_status", "std::vector<double>", sim_GenParticle_status);
 		OutputTree->Branch("GenParticle_time", "std::vector<double>", sim_GenParticle_time);
 		OutputTree->Branch("GenParticle_x", "std::vector<double>", sim_GenParticle_x);
 		OutputTree->Branch("GenParticle_y", "std::vector<double>", sim_GenParticle_y);
 		OutputTree->Branch("GenParticle_z", "std::vector<double>", sim_GenParticle_z);
 		OutputTree->Branch("GenParticle_energy", "std::vector<double>", sim_GenParticle_energy);
 		OutputTree->Branch("GenParticle_px", "std::vector<double>", sim_GenParticle_px);
 		OutputTree->Branch("GenParticle_py", "std::vector<double>", sim_GenParticle_py);
 		OutputTree->Branch("GenParticle_pz", "std::vector<double>", sim_GenParticle_pz);
//		OutputTree->Branch("GenParticle_p_mag", &p);
 	//	OutputTree->Branch("GenParticle_mo1", "std::vector<double>", sim_GenParticle_mo1);
 	//	OutputTree->Branch("GenParticle_mo2", "std::vector<double>", sim_GenParticle_mo2);
 	//	OutputTree->Branch("GenParticle_dau1", "std::vector<double>", sim_GenParticle_dau1);
 	//	OutputTree->Branch("GenParticle_dau2", "std::vector<double>", sim_GenParticle_dau2);
 		OutputTree->Branch("GenParticle_mass", "std::vector<double>", sim_GenParticle_mass);
 //		OutputTree->Branch("GenParticle_pt", "std::vector<double>", sim_GenParticle_pt);
 //		OutputTree->Branch("GenParticle_eta", "std::vector<double>", sim_GenParticle_eta);
 //		OutputTree->Branch("GenParticle_phi", "std::vector<double>", sim_GenParticle_phi);
 //		OutputTree->Branch("COSMIC_EVENT_ID", "std::vector<double>", sim_COSMIC_EVENT_ID);
 //		OutputTree->Branch("COSMIC_CORE_X", "std::vector<double>", sim_COSMIC_CORE_X);
 //		OutputTree->Branch("COSMIC_CORE_Y", "std::vector<double>", sim_COSMIC_CORE_Y);
 //		OutputTree->Branch("COSMIC_GEN_PRIMARY_ENERGY", "std::vector<double>", sim_COSMIC_GEN_PRIMARY_ENERGY);
 //		OutputTree->Branch("COSMIC_GEN_THETA", "std::vector<double>", sim_COSMIC_GEN_THETA);
 //		OutputTree->Branch("COSMIC_GEN_PHI", "std::vector<double>", sim_COSMIC_GEN_PHI);
 //		OutputTree->Branch("COSMIC_GEN_FIRST_HEIGHT", "std::vector<double>", sim_COSMIC_GEN_FIRST_HEIGHT);
 //		OutputTree->Branch("COSMIC_GEN_ELECTRON_COUNT", "std::vector<double>", sim_COSMIC_GEN_ELECTRON_COUNT);
 //		OutputTree->Branch("COSMIC_GEN_MUON_COUNT", "std::vector<double>", sim_COSMIC_GEN_MUON_COUNT);
 //		OutputTree->Branch("COSMIC_GEN_HADRON_COUNT", "std::vector<double>", sim_COSMIC_GEN_HADRON_COUNT);
 //		OutputTree->Branch("COSMIC_GEN_PRIMARY_ID", "std::vector<double>", sim_COSMIC_GEN_PRIMARY_ID);
 		OutputTree->Branch("EXTRA_11", "std::vector<double>", sim_EXTRA_11);
 		OutputTree->Branch("EXTRA_12", "std::vector<double>", sim_EXTRA_12);
 		OutputTree->Branch("EXTRA_13", "std::vector<double>", sim_EXTRA_13);
 		OutputTree->Branch("EXTRA_14", "std::vector<double>", sim_EXTRA_14);
 		OutputTree->Branch("EXTRA_15", "std::vector<double>", sim_EXTRA_15);

 		OutputTree->Branch("NumVertices", &numvertices, "NumVertices/D");
      	OutputTree->Branch("Vertex_numTracks", &vertex_numTracks);
        OutputTree->Branch("Vertex_cosOpeningAngle", &vertex_cosOpeningAngle);
      	OutputTree->Branch("Vertex_t", &vertex_t);
      	OutputTree->Branch("Vertex_x", &vertex_x);
      	OutputTree->Branch("Vertex_y", &vertex_y);
      	OutputTree->Branch("Vertex_z", &vertex_z);
      	OutputTree->Branch("Vertex_ErrorT", &vertex_t_error);
      	OutputTree->Branch("Vertex_ErrorX", &vertex_x_error);
      	OutputTree->Branch("Vertex_ErrorY", &vertex_y_error);
      	OutputTree->Branch("Vertex_ErrorZ", &vertex_z_error);
      	OutputTree->Branch("Vertex_chi2", &vertex_chi2);
      	OutputTree->Branch("Vertex_chi2PerNdof", &vertex_chi2_per_dof);
      	OutputTree->Branch("Vertex_chi2PValue", &vertex_chi2_p_value);
        OutputTree->Branch("Vertex_trackIndices", &vertex_track_indices);

//        OutputTree->Branch("Vertex_k_t", &vertex_k_t);
//        OutputTree->Branch("Vertex_k_x", &vertex_k_x);
//        OutputTree->Branch("Vertex_k_y", &vertex_k_y);
//        OutputTree->Branch("Vertex_k_z", &vertex_k_z);
//        OutputTree->Branch("Vertex_k_trackIndices", &vertex_track_k_indices);

	OutputTree->Branch("Vertex_k_m_t", &vertex_k_m_t);
        OutputTree->Branch("Vertex_k_m_x", &vertex_k_m_x);
        OutputTree->Branch("Vertex_k_m_y", &vertex_k_m_y);
        OutputTree->Branch("Vertex_k_m_z", &vertex_k_m_z);
      	OutputTree->Branch("Vertex_k_m_ErrorT", &vertex_t_k_m_error);
      	OutputTree->Branch("Vertex_k_m_ErrorX", &vertex_x_k_m_error);
      	OutputTree->Branch("Vertex_k_m_ErrorY", &vertex_y_k_m_error);
      	OutputTree->Branch("Vertex_k_m_ErrorZ", &vertex_z_k_m_error);
        OutputTree->Branch("Vertex_k_m_trackIndices", &vertex_track_k_m_indices);
	OutputTree->Branch("NumVertices_k_m", &numvertices_k_m, "NumVertices/D");

      	OutputTree->Branch("NumTracks", &numtracks, "NumTracks/D");
      	OutputTree->Branch("Track_numHits", &track_numHits);
      	OutputTree->Branch("Track_t0", &track_t);
      	OutputTree->Branch("Track_x0", &track_x);
      	OutputTree->Branch("Track_y0", &track_y);
      	OutputTree->Branch("Track_z0", &track_z);
      	OutputTree->Branch("Track_velX", &track_vx);
      	OutputTree->Branch("Track_velY", &track_vy);
      	OutputTree->Branch("Track_velZ", &track_vz);
      	OutputTree->Branch("Track_ErrorT0", &track_t_error);
      	OutputTree->Branch("Track_ErrorX0", &track_x_error);
      	OutputTree->Branch("Track_ErrorY0", &track_y_error);
      	OutputTree->Branch("Track_ErrorZ0", &track_z_error);
      	OutputTree->Branch("Track_ErrorVx", &track_vx_error);
      	OutputTree->Branch("Track_ErrorVy", &track_vy_error);
      	OutputTree->Branch("Track_ErrorVz", &track_vz_error);
      	OutputTree->Branch("Track_chi2", &track_chi2);
      	OutputTree->Branch("Track_chi2PerNdof", &track_chi2_per_dof);
      	OutputTree->Branch("Track_chi2PValue", &track_chi2_p_value);
      	OutputTree->Branch("Track_beta", &track_beta);
      	OutputTree->Branch("Track_beta_err", &track_beta_error);
      	OutputTree->Branch("Track_ErrorBeta", &track_beta_error);
      	OutputTree->Branch("Track_angle", &track_angle);
       	OutputTree->Branch("Track_ErrorAngle", &track_angle_error);
       	OutputTree->Branch("Track_detCount", &unique_detector_count);
      	OutputTree->Branch("Track_expectedHitLayer", &track_expected_hit_layer);
      	OutputTree->Branch("Track_missingHitLayer", &track_missing_hit_layer);
        OutputTree->Branch("Track_hitIndices", &track_hit_indices);
        OutputTree->Branch("track_ipDistance", &track_distance_to_ip);

	OutputTree->Branch("local_chi_f", &chi_f_first);
        OutputTree->Branch("local_chi_s", &chi_s_first);
	OutputTree->Branch("Track_k_filterchi", &track_filter_chi);
	OutputTree->Branch("Track_k_smoothchi", &track_smooth_chi);
	OutputTree->Branch("Track_k_numHits", &track_k_numHits);
	OutputTree->Branch("Track_k_chi2PerNdof", &track_k_chi2_per_dof);
	OutputTree->Branch("Track_k_smooth_chi_sum", &track_k_smooth_chi_sum);
/*
	OutputTree->Branch("Track_k_hitIndices", &track_hit_k_indices);
	OutputTree->Branch("Track_k_velX", &track_k_vx);
	OutputTree->Branch("Track_k_velY", &track_k_vy);
	OutputTree->Branch("Track_k_velZ", &track_k_vz);
	OutputTree->Branch("Track_k_beta", &track_k_beta);
	OutputTree->Branch("Track_k_beta_err", &track_k_beta_err);
	OutputTree->Branch("Track_k_x0", &track_k_x);
	OutputTree->Branch("Track_k_y0", &track_k_y);
	OutputTree->Branch("Track_k_z0", &track_k_z);
	OutputTree->Branch("Track_k_t0", &track_k_t);
 	OutputTree->Branch("NumTracks_k", &NumTracks_k);
	OutputTree->Branch("x_estimates", &x_estimates);
	OutputTree->Branch("y_estimates", &y_estimates);
	OutputTree->Branch("z_estimates", &z_estimates);
*/
	OutputTree->Branch("x_std_scat_per_m", &x_scat);
	OutputTree->Branch("z_std_scat_per_m", &z_scat);
	OutputTree->Branch("track_pdgs", &track_pdgs);
	OutputTree->Branch("track_ids", &track_ids);
	OutputTree->Branch("track_k_opening_angle", &track_k_openingangle);

	OutputTree->Branch("Track_k_m_velX", &track_k_m_vx);
	OutputTree->Branch("Track_k_m_velY", &track_k_m_vy);
	OutputTree->Branch("Track_k_m_velZ", &track_k_m_vz);
	OutputTree->Branch("Track_k_m_x0", &track_k_m_x);
	OutputTree->Branch("Track_k_m_y0", &track_k_m_y);
	OutputTree->Branch("Track_k_m_z0", &track_k_m_z);
	OutputTree->Branch("Track_k_m_t0", &track_k_m_t);
	OutputTree->Branch("NumTracks_k_m", &NumTracks_k_m);
	OutputTree->Branch("Track_k_m_hitIndices", &track_hit_k_m_indices);
	OutputTree->Branch("Track_k_m_expected_hit_layer", &track_k_m_expected_hit_layer);
	OutputTree->Branch("x_estimates_m", &x_estimates_m);
	OutputTree->Branch("y_estimates_m", &y_estimates_m);
	OutputTree->Branch("z_estimates_m", &z_estimates_m);
/*
	OutputTree->Branch("vertex_vx_m", &q_s_x_m);
	OutputTree->Branch("vertex_vy_m", &q_s_y_m);
	OutputTree->Branch("vertex_vz_m", &q_s_z_m);
	OutputTree->Branch("vertex_k_f_beta", &vertex_k_f_beta);
	OutputTree->Branch("vertex_k_s_beta", &vertex_k_s_beta);
	OutputTree->Branch("vertex_k_s_beta_err", &vertex_k_s_beta_err);
	OutputTree->Branch("king_move_inds", &king_move_inds);
*/


	}


  //____________________________________________________________________________________________

 //___Make Sim Branches_________________________________________________________________________
 	Double_t sim_numhits;
 	std::vector<double> *sim_hit_e = nullptr;
 	std::vector<double> *sim_hit_t = nullptr;
 	std::vector<double> *sim_hit_detId = nullptr;
 	std::vector<double> *sim_hit_particlePdgId = nullptr;
 	std::vector<double> *sim_hit_G4TrackId = nullptr;
 	std::vector<double> *sim_hit_G4ParentTrackId = nullptr;
 	std::vector<double> *sim_hit_x = nullptr;
 	std::vector<double> *sim_hit_y = nullptr;
 	std::vector<double> *sim_hit_z = nullptr;
 	std::vector<double> *sim_hit_particleEnergy = nullptr;
 	std::vector<double> *sim_hit_px = nullptr;
 	std::vector<double> *sim_hit_py = nullptr;
 	std::vector<double> *sim_hit_pz = nullptr;
 	std::vector<double> *sim_hit_weight = nullptr;
 	Double_t sim_NumGenParticles;
 	std::vector<double> *sim_GenParticle_index = nullptr;
 	std::vector<double> *sim_GenParticle_G4index = nullptr;
 	std::vector<double> *sim_GenParticle_pdgid = nullptr;
 	std::vector<double> *sim_GenParticle_status = nullptr;
 	std::vector<double> *sim_GenParticle_time = nullptr;
 	std::vector<double> *sim_GenParticle_x = nullptr;
 	std::vector<double> *sim_GenParticle_y = nullptr;
 	std::vector<double> *sim_GenParticle_z = nullptr;
 	std::vector<double> *sim_GenParticle_energy = nullptr;
 	std::vector<double> *sim_GenParticle_px = nullptr;
 	std::vector<double> *sim_GenParticle_py = nullptr;
 	std::vector<double> *sim_GenParticle_pz = nullptr;
 	std::vector<double> *sim_GenParticle_mo1 = nullptr;
 	std::vector<double> *sim_GenParticle_mo2 = nullptr;
 	std::vector<double> *sim_GenParticle_dau1 = nullptr;
 	std::vector<double> *sim_GenParticle_dau2 = nullptr;
 	std::vector<double> *sim_GenParticle_mass = nullptr;
 	std::vector<double> *sim_GenParticle_pt = nullptr;
 	std::vector<double> *sim_GenParticle_eta = nullptr;
 	std::vector<double> *sim_GenParticle_phi = nullptr;
	std::vector<double> *sim_COSMIC_EVENT_ID = nullptr;
 	std::vector<double> *sim_COSMIC_CORE_X = nullptr;
 	std::vector<double> *sim_COSMIC_CORE_Y = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_PRIMARY_ENERGY = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_THETA = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_PHI = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_FIRST_HEIGHT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_ELECTRON_COUNT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_MUON_COUNT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_HADRON_COUNT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_PRIMARY_ID = nullptr;
 	std::vector<double> *sim_EXTRA_11 = nullptr;
 	std::vector<double> *sim_EXTRA_12 = nullptr;
 	std::vector<double> *sim_EXTRA_13 = nullptr;
 	std::vector<double> *sim_EXTRA_14 = nullptr;
 	std::vector<double> *sim_EXTRA_15 = nullptr;

 		//__Make Vertex Branches________________________________________________________________________
  	std::vector<int> vertex_numTracks;
    std::vector<double> vertex_cosOpeningAngle;
  	std::vector<double> vertex_chi2;
  	std::vector<double> vertex_chi2_per_dof;
 	std::vector<double> vertex_chi2_p_value;
  	std::vector<double> vertex_t;
  	std::vector<double> vertex_x;
  	std::vector<double> vertex_y;
  	std::vector<double> vertex_z;
  	std::vector<double> vertex_t_error;
  	std::vector<double> vertex_x_error;
  	std::vector<double> vertex_y_error;
  	std::vector<double> vertex_z_error;
    std::vector<int> vertex_track_indices;
  	Double_t numvertices;

        std::vector<double> vertex_k_t;
        std::vector<double> vertex_k_x;
        std::vector<double> vertex_k_y;
        std::vector<double> vertex_k_z;
	std::vector<int> vertex_track_k_indices;

	std::vector<double> q_s_x_m;
	std::vector<double> q_s_y_m;
	std::vector<double> q_s_z_m;
	std::vector<double> vertex_k_f_beta;
	std::vector<double> vertex_k_s_beta;
	std::vector<double> vertex_k_s_beta_err;

std::vector<double>	vertex_k_m_t;
std::vector<double>	vertex_k_m_x;
std::vector<double>	vertex_k_m_y;
std::vector<double>	vertex_k_m_z;
  	std::vector<double> vertex_t_k_m_error;
  	std::vector<double> vertex_x_k_m_error;
  	std::vector<double> vertex_y_k_m_error;
  	std::vector<double> vertex_z_k_m_error;
std::vector<double>	vertex_track_k_m_indices;
Double_t numvertices_k_m;

  //__Make Track Branches_________________________________________________________________________
  	std::vector<double> track_numHits;
  	std::vector<double> track_chi2;
  	std::vector<double> track_chi2_per_dof;
  	std::vector<double> track_chi2_p_value;
  	std::vector<double> track_beta;
  	std::vector<double> track_beta_error;
  	std::vector<double> track_angle;
  	std::vector<double> track_angle_error;
  	std::vector<double> unique_detector_count;
  	std::vector<double> track_t;
  	std::vector<double> track_x;
  	std::vector<double> track_y;
  	std::vector<double> track_z;
  	std::vector<double> track_vx;
  	std::vector<double> track_vy;
  	std::vector<double> track_vz;
  	std::vector<double> track_t_error;
  	std::vector<double> track_x_error;
  	std::vector<double> track_y_error;
  	std::vector<double> track_z_error;
  	std::vector<double> track_vx_error;
  	std::vector<double> track_vy_error;
  	std::vector<double> track_vz_error;
  	std::vector<int> track_expected_hit_layer;
  	std::vector<double> track_missing_hit_layer;
        std::vector<int> track_hit_indices;
        std::vector<double> track_distance_to_ip;
        std::vector<double> track_filter_chi;
        std::vector<double> track_smooth_chi;
  	Double_t numtracks;
		std::vector<double> track_k_numHits;
		std::vector<double> track_k_chi2_per_dof;
		std::vector<double> track_k_smooth_chi_sum;
		std::vector<int> track_hit_k_indices;

	std::vector<double> chi_f_first;
	std::vector<double> chi_s_first;
	std::vector<double> track_k_vx;
	std::vector<double> track_k_vy;
	std::vector<double> track_k_vz;
	std::vector<double> track_k_beta;
	std::vector<double> track_k_beta_err;
	std::vector<double> track_k_t;
	std::vector<double> track_k_x;
	std::vector<double> track_k_y;
	std::vector<double> track_k_z;
	int NumTracks_k;
	std::vector<double> x_estimates;
        std::vector<double> y_estimates;
        std::vector<double> z_estimates;
	std::vector<double> x_scat;
	std::vector<double> z_scat;
	std::vector<int> track_pdgs;
	std::vector<int> track_ids;
	std::vector<double> track_k_openingangle;


	std::vector<double> track_k_m_vx;
	std::vector<double> track_k_m_vy;
	std::vector<double> track_k_m_vz;
	std::vector<double> track_k_m_t;
	std::vector<double> track_k_m_x;
	std::vector<double> track_k_m_y;
	std::vector<double> track_k_m_z;
	int NumTracks_k_m;
	std::vector<int> track_hit_k_m_indices;
	std::vector<double> x_estimates_m;
        std::vector<double> y_estimates_m;
        std::vector<double> z_estimates_m;
	std::vector<double> track_k_m_expected_hit_layer;
	std::vector<int> king_move_inds;

  //___Make Digi Branches_____________________________________________________________________
  	std::vector<double> digi_hit_t;
  	std::vector<double> digi_hit_x;
  	std::vector<double> digi_hit_y;
  	std::vector<double> digi_hit_z;
  	std::vector<double> digi_hit_e;
  	std::vector<double> digi_hit_px;
  	std::vector<double> digi_hit_py;
  	std::vector<double> digi_hit_pz;
    std::vector<double> digi_particle_energy;
    std::vector<int> digi_pdg;
  	std::vector<int> digi_hit_indices;
  	std::vector<int> Digi_numHits;

}; //class TreeHandler

template<class digi_hit>
void TreeHandler::ExportDigis(std::vector<digi_hit*> digi_list){
      digi_hit_indices.clear();
      digi_hit_t.clear();
      digi_hit_x.clear();
      digi_hit_y.clear();
      digi_hit_z.clear();
      digi_hit_e.clear();
      digi_hit_px.clear();
      digi_hit_py.clear();
      digi_hit_pz.clear();
      Digi_numHits.clear();
      digi_particle_energy.clear();
      digi_pdg.clear();

      for (auto digi : digi_list){
        Digi_numHits.push_back(digi->hits.size());
        digi_hit_t.push_back(digi->t);
        digi_hit_x.push_back(digi->x);
        digi_hit_y.push_back(digi->y);
        digi_hit_z.push_back(digi->z);
        digi_hit_e.push_back(digi->e);
        digi_hit_px.push_back(digi->px);
        digi_hit_py.push_back(digi->py);
        digi_hit_pz.push_back(digi->pz);
        digi_particle_energy.push_back(digi->particle_energy);
        digi_pdg.push_back(digi->pdg);
        for (auto hit : digi->hits){
          digi_hit_indices.push_back(hit->index);
        }
      }

 }

template<class Track>
void TreeHandler::ExportTracks_k(std::vector<Track*> track_list){

    //track_k_numhits.clear();
    track_filter_chi.clear();
    track_smooth_chi.clear();
    track_k_numHits.clear();
    track_k_chi2_per_dof.clear();
    track_k_smooth_chi_sum.clear();
    track_hit_k_indices.clear();
    track_k_vx.clear();
    track_k_vy.clear();
    track_k_vz.clear();
    track_k_x.clear();
    track_k_y.clear();
    track_k_z.clear();
    track_k_t.clear();
    x_estimates.clear();
    y_estimates.clear();
    z_estimates.clear();
    x_scat.clear();
    z_scat.clear();
    track_pdgs.clear();
    track_ids.clear();
    track_k_openingangle.clear();
    track_k_beta.clear();
    track_k_beta_err.clear();

    NumTracks_k = track_list.size(); // this is number of tracks not hits!

//find two tracks
	// 	for (auto tr1 : track_list){
	// 		for (auto tr2 : track_list){
	//
	// 	// if (NumTracks_k == 2){
	// 	// 	track_k_openingangle.push_back(track_list[0]->direction()^track_list[1]->direction());
	// 	// }
	// }
	// }
	for (int i = 0; i < NumTracks_k; i++) {
		for (int j = i+1; j < NumTracks_k; j++) {
		auto tr1 =  track_list[i];
		auto tr2 =  track_list[j];
		track_k_openingangle.push_back(tr1->direction()^tr2->direction());
	}
}

    for (auto tr : track_list){
			track_k_numHits.push_back( (tr->hits).size() );
			track_k_chi2_per_dof.push_back(tr->chi2_per_dof());
      for (auto chi : tr->chi_f) track_filter_chi.push_back( static_cast<double>(chi) );
      track_filter_chi.push_back(-1);

//        std::cout << " Tree Chis are: " << std::endl;
      for (auto chi : tr->chi_s) {track_smooth_chi.push_back( static_cast<double>(chi) );}
  //      std::cout << ", " << chi;}
//	std::cout << std::endl;

      track_smooth_chi.push_back(-1);

			double chi_sum = 0;
			for (auto chi : tr->chi_s) {
				chi_sum += chi;
			}
			chi_sum = chi_sum / (4.0*tr->chi_s.size() - 6.0);
			track_k_smooth_chi_sum.push_back(chi_sum);

			for (auto hit : tr->hits) { track_hit_k_indices.push_back(hit->index); }
			track_hit_k_indices.push_back(-1.);

			for (auto est : tr->estimate_list) {
				x_estimates.push_back(est[0]);
				y_estimates.push_back(est[1]);
				z_estimates.push_back(est[2]);
			}
			x_estimates.push_back(-1.0);
			y_estimates.push_back(-1.0);
			z_estimates.push_back(-1.0);

			for (auto std : tr->x_scats) x_scat.push_back(std);
			for (auto std : tr->z_scats) z_scat.push_back(std);
			//x_scat.push_back(-1.0);
			//z_scat.push_back(-1.0);

			track_k_vx.push_back(tr->vx);
 			track_k_vy.push_back(tr->vy);
			track_k_vz.push_back(tr->vz);

			double v = std::sqrt(std::pow(tr->vx,2) + std::pow(tr->vy,2) + std::pow(tr->vz,2));
			track_k_beta.push_back(v / constants::c);

			Eigen::MatrixXd R = tr->P_s;
			Eigen::MatrixXd D(3,3);
			D << R(3,3), R(3,4), R(3,5),
		             R(4,3), R(4,4), R(4,5),
       			     R(5,3), R(5,4), R(5,5);

			Eigen::VectorXd q(3);
			q << tr->vx, tr->vy, tr->vz;

			track_k_beta_err.push_back((q.transpose() * D * q)(0) / (v * v));

			track_k_t.push_back(tr->t0);
			track_k_x.push_back(tr->x0);
			track_k_y.push_back(tr->y0);
			track_k_z.push_back(tr->z0);

			for (auto hit : tr->hits) {
				track_ids.push_back(hit->min_track_id);
				track_pdgs.push_back(hit->pdg);
			}
			track_ids.push_back(-1);
			track_pdgs.push_back(-1);

    }
}

template<class Track>
void TreeHandler::ExportTracks_k_m(std::vector<Track*> track_list){
	track_k_m_vx.clear();
	track_k_m_vy.clear();
	track_k_m_vz.clear();
	track_k_m_x.clear();
	track_k_m_y.clear();
	track_k_m_z.clear();
	track_k_m_t.clear();
	track_hit_k_m_indices.clear();

	x_estimates_m.clear();
	y_estimates_m.clear();
	z_estimates_m.clear();
	track_k_m_expected_hit_layer.clear();
	king_move_inds.clear();
	NumTracks_k_m = track_list.size(); // this is number of tracks not hits!


    for (auto tr : track_list){
	track_k_m_vx.push_back(tr->vx);
	track_k_m_vy.push_back(tr->vy);
	track_k_m_vz.push_back(tr->vz);
	track_k_m_t.push_back(tr->t0);
	track_k_m_x.push_back(tr->x0);
	track_k_m_y.push_back(tr->y0);
	track_k_m_z.push_back(tr->z0);

	//if (tr->king_move_inds.size() != 0) std::cout << "\n Treehandler ";
	//std::cout << "kingly Tree ";
	for (auto ind : tr->king_move_inds) {
		king_move_inds.push_back(ind);
		//std::cout << ind << ", ";
	}
	if (tr->king_move_inds.size() != 0) king_move_inds.push_back(-1);
	//std::cout << std::endl;

	for (auto hit : tr->hits) {
		track_hit_k_m_indices.push_back(hit->index);
	}
	track_hit_k_m_indices.push_back(-1.);

	for (auto exp_layer : tr->expected_layers) {
		track_k_m_expected_hit_layer.push_back( static_cast<double>(exp_layer) );
	}
	track_k_m_expected_hit_layer.push_back(-1);

	//for (auto q_s : tr->q_s) {
	//}



	for (auto est : tr->estimate_list) {
		x_estimates_m.push_back(est[0]);
		y_estimates_m.push_back(est[1]);
		z_estimates_m.push_back(est[2]);
	}
	x_estimates_m.push_back(-1.0);
	y_estimates_m.push_back(-1.0);
	z_estimates_m.push_back(-1.0);


}

}

void TreeHandler::Export_kalman_info(std::vector<double> local_chi_f, std::vector<double> local_chi_s) {

	chi_f_first.clear();
	chi_s_first.clear();

	for (auto chi : local_chi_f) chi_f_first.push_back(chi);
	for (auto chi : local_chi_s) chi_s_first.push_back(chi);
}


template<class Track>
void TreeHandler::ExportTracks(std::vector<Track*> track_list){

    track_numHits.clear();
    track_chi2.clear();
    track_chi2_per_dof.clear();
    track_chi2_p_value.clear();
    track_beta.clear();
    track_beta_error.clear();
    track_angle.clear();
    track_angle_error.clear();
    unique_detector_count.clear();
    track_t.clear();
    track_x.clear();
    track_y.clear();
    track_z.clear();
    track_vx.clear();
    track_vy.clear();
    track_vz.clear();
    track_t_error.clear();
    track_x_error.clear();
    track_y_error.clear();
    track_z_error.clear();
    track_vx_error.clear();
    track_vy_error.clear();
    track_vz_error.clear();
    track_hit_indices.clear();
    track_missing_hit_layer.clear();
    track_distance_to_ip.clear();
    track_expected_hit_layer.clear();

    numtracks = track_list.size();
    for (auto tr : track_list){
      track_numHits.push_back( (tr->hits).size() );
      track_chi2.push_back(tr->chi2());
      track_chi2_per_dof.push_back(tr->chi2_per_dof());
      track_beta.push_back(tr->beta());
      track_beta_error.push_back(tr->beta_err());
      track_t.push_back(tr->t0);
      track_x.push_back(tr->x0);
      track_y.push_back(tr->y0);
      track_z.push_back(tr->z0);
      track_vx.push_back(tr->vx);
      track_vy.push_back(tr->vy);
      track_vz.push_back(tr->vz);
      track_x_error.push_back(tr->ex0);
      track_y_error.push_back(tr->ey0);
      track_z_error.push_back(tr->ez0);
      track_vx_error.push_back(tr->evx);
      track_vy_error.push_back(tr->evy);
      track_vz_error.push_back(tr->evz);
      track_angle.push_back(tr->cos_angle_from_ip());
      track_distance_to_ip.push_back(tr->shortDistance());

      for (auto missing_layer : tr->_missing_layers) track_missing_hit_layer.push_back( static_cast<double>(missing_layer) );
      track_missing_hit_layer.push_back(-1);

      for (auto exp_layer : tr->expected_layers) track_expected_hit_layer.push_back( static_cast<double>(exp_layer) );
      track_expected_hit_layer.push_back(-1);

      for (auto hit : tr->hits) { track_hit_indices.push_back(hit->index); }
      track_hit_indices.push_back(-1.);

    }

}

template<typename vertex>
void TreeHandler::ExportVertices_k(std::vector<vertex*> vertices){

	vertex_k_t.clear();
	vertex_k_x.clear();
	vertex_k_y.clear();
	vertex_k_z.clear();
	vertex_track_k_indices.clear();

	for (auto v : vertices) {

		vertex_k_t.push_back(v->t);
		vertex_k_x.push_back(v->x);
		vertex_k_y.push_back(v->y);
		vertex_k_z.push_back(v->z);

		for (auto tr_index : v->track_indices){
        	vertex_track_k_indices.push_back(tr_index);
        	}

	        vertex_track_k_indices.push_back(-1);
	}

}

template<typename vertex>
void TreeHandler::ExportVertices_k_m(std::vector<vertex*> vertices){

	vertex_k_m_t.clear();
	vertex_k_m_x.clear();
	vertex_k_m_y.clear();
	vertex_k_m_z.clear();
	vertex_track_k_m_indices.clear();

	q_s_x_m.clear();
	q_s_y_m.clear();
	q_s_z_m.clear();
	vertex_k_f_beta.clear();
	vertex_k_s_beta.clear();
	vertex_k_s_beta_err.clear();

	vertex_t_k_m_error.clear();
	vertex_x_k_m_error.clear();
	vertex_y_k_m_error.clear();
	vertex_z_k_m_error.clear();
	numvertices_k_m = vertices.size();

	//std::cout << "export test 1 " << std::endl;

	for (auto v : vertices) {

		//std::cout << "export test 2 " << std::endl;

		vertex_k_m_t.push_back(v->t);
		vertex_k_m_x.push_back(v->x);
		vertex_k_m_y.push_back(v->y);
		vertex_k_m_z.push_back(v->z);

		//std::cout << "export test 3 " << std::endl;
		vertex_x_k_m_error.push_back(sqrt(v->CovMatrix()[0][0]));
		vertex_y_k_m_error.push_back(sqrt(v->CovMatrix()[1][1]));
		vertex_z_k_m_error.push_back(sqrt(v->CovMatrix()[2][2]));
		vertex_t_k_m_error.push_back(sqrt(v->CovMatrix()[3][3]));

		for (auto q_f : v->q_f) vertex_k_f_beta.push_back(q_f.norm() / constants::c);
		vertex_k_f_beta.push_back(-1.0);

		for (int i=0; i < v->q_s.size(); i++) {
			q_s_x_m.push_back(v->q_s[i][0]);
			q_s_y_m.push_back(v->q_s[i][1]);
			q_s_z_m.push_back(v->q_s[i][2]);

			//double pull = (v->q_s[i].norm() - constants::c) * v->q_s[i].squaredNorm() /
			//	(v->q_s[i].transpose() * v->D_s[i] * v->q_s[i])(0);
			//vertex_k_beta.push_back(pull);

			double vertex_k_beta_err = (v->q_s[i].transpose() * v->D_s[i] * v->q_s[i])(0) / v->q_s[i].squaredNorm();

			vertex_k_s_beta_err.push_back(vertex_k_beta_err);

			vertex_k_s_beta.push_back((v->q_s[i]).norm() / constants::c);
		}
		q_s_x_m.push_back(-1.0);
		q_s_y_m.push_back(-1.0);
		q_s_z_m.push_back(-1.0);
		vertex_k_s_beta.push_back(-1.0);

		for (auto tr_index : v->track_indices){
        		vertex_track_k_m_indices.push_back(tr_index);
        	}

		//std::cout << "export test 4 " << std::endl;

	        vertex_track_k_m_indices.push_back(-1);
	}

}

template<typename vertex>
void TreeHandler::ExportVertices(std::vector<vertex*> vertices){

    vertex_numTracks.clear();
    vertex_cosOpeningAngle.clear();
    vertex_chi2.clear();
    vertex_chi2_per_dof.clear();
    vertex_chi2_p_value.clear();
    vertex_t.clear();
    vertex_x.clear();
    vertex_y.clear();
    vertex_z.clear();
    vertex_t_error.clear();
    vertex_x_error.clear();
    vertex_y_error.clear();
    vertex_z_error.clear();
    vertex_track_indices.clear();


    numvertices = vertices.size();


    for (auto v : vertices){
      vertex_numTracks.push_back(v->track_indices.size());
      vertex_cosOpeningAngle.push_back(v->cos_opening_angle);
      vertex_t.push_back(v->t);
      vertex_x.push_back(v->x);
      vertex_y.push_back(v->y);
      vertex_z.push_back(v->z);

      vertex_x_error.push_back(sqrt(v->CovMatrix()[0][0]));
      vertex_y_error.push_back(sqrt(v->CovMatrix()[1][1]));
      vertex_z_error.push_back(sqrt(v->CovMatrix()[2][2]));
      vertex_t_error.push_back(sqrt(v->CovMatrix()[3][3]));

      vertex_chi2_per_dof.push_back(v->merit());


      for (auto tr_index : v->track_indices){
        vertex_track_indices.push_back(tr_index);
      }

      vertex_track_indices.push_back(-1);

    }

}

#endif
