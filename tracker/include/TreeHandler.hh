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

  //void Export_kalman_info(std::vector<double>, std::vector<double>);

  template<class digi_hit>
  void ExportDigis(std::vector<digi_hit*>, long long int digi_seed);

  //template<class track>
  //void ExportTracks_k(std::vector<track*>);

	template<class track>
  void ExportTracks_k_m(std::vector<track*>);

	template<class track>
  void ExportTracks_c_b(std::vector<track*>);

  	template<class track>
  void ExportTracks(std::vector<track*>);

  	template<typename vertex>
  void ExportVertices(std::vector<vertex*>);

	 template<typename vertex>
  void ExportVertices_c_b(std::vector<vertex*>);

  //template<typename vertex>
  //void ExportVertices_k(std::vector<vertex*>);

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
        OutputTree->Branch("Digi_seed", &digi_seed, "Digi_seed/L");
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

	OutputTree->Branch("Vertex_k_m_t", &vertex_k_m_t);
        OutputTree->Branch("Vertex_k_m_x", &vertex_k_m_x);
        OutputTree->Branch("Vertex_k_m_y", &vertex_k_m_y);
        OutputTree->Branch("Vertex_k_m_z", &vertex_k_m_z);
      	OutputTree->Branch("Vertex_k_m_ErrorT", &vertex_t_k_m_error);
      	OutputTree->Branch("Vertex_k_m_ErrorX", &vertex_x_k_m_error);
      	OutputTree->Branch("Vertex_k_m_ErrorY", &vertex_y_k_m_error);
      	OutputTree->Branch("Vertex_k_m_ErrorZ", &vertex_z_k_m_error);
	OutputTree->Branch("vertex_k_m_chi2", &vertex_k_m_chi2_per_dof);
        OutputTree->Branch("Vertex_k_m_trackIndices", &vertex_track_k_m_indices);
	OutputTree->Branch("NumVertices_k_m", &numvertices_k_m, "NumVertices/D");

	OutputTree->Branch("Track_k_m_velX", &track_k_m_vx);
	OutputTree->Branch("Track_k_m_velY", &track_k_m_vy);
	OutputTree->Branch("Track_k_m_velZ", &track_k_m_vz);
	OutputTree->Branch("Track_k_m_x0", &track_k_m_x);
	OutputTree->Branch("Track_k_m_y0", &track_k_m_y);
	OutputTree->Branch("Track_k_m_z0", &track_k_m_z);
	OutputTree->Branch("Track_k_m_t0", &track_k_m_t);

      	OutputTree->Branch("Track_k_m_ErrorT0", &track_k_m_t_error);
      	OutputTree->Branch("Track_k_m_ErrorX0", &track_k_m_x_error);
      	OutputTree->Branch("Track_k_m_ErrorY0", &track_k_m_y_error);
      	OutputTree->Branch("Track_k_m_ErrorZ0", &track_k_m_z_error);
      	OutputTree->Branch("Track_k_m_ErrorVx", &track_k_m_vx_error);
      	OutputTree->Branch("Track_k_m_ErrorVy", &track_k_m_vy_error);
      	OutputTree->Branch("Track_k_m_ErrorVz", &track_k_m_vz_error);

	OutputTree->Branch("Track_k_m_cov_x_t", &track_k_m_cov_x_t);
	OutputTree->Branch("Track_k_m_cov_x_z", &track_k_m_cov_x_z);
	OutputTree->Branch("Track_k_m_cov_x_vx", &track_k_m_cov_x_vx);
	OutputTree->Branch("Track_k_m_cov_x_vy", &track_k_m_cov_x_vy);
	OutputTree->Branch("Track_k_m_cov_x_vz", &track_k_m_cov_x_vz);

	OutputTree->Branch("Track_k_m_cov_t_z", &track_k_m_cov_t_z);
	OutputTree->Branch("Track_k_m_cov_t_vx", &track_k_m_cov_t_vx);
	OutputTree->Branch("Track_k_m_cov_t_vy", &track_k_m_cov_t_vy);
	OutputTree->Branch("Track_k_m_cov_t_vz", &track_k_m_cov_t_vz);

	OutputTree->Branch("Track_k_m_cov_z_vx", &track_k_m_cov_z_vx);
	OutputTree->Branch("Track_k_m_cov_z_vy", &track_k_m_cov_z_vy);
	OutputTree->Branch("Track_k_m_cov_z_vz", &track_k_m_cov_z_vz);

	OutputTree->Branch("Track_k_m_cov_vx_vy", &track_k_m_cov_vx_vy);
	OutputTree->Branch("Track_k_m_cov_vx_vz", &track_k_m_cov_vx_vz);

	OutputTree->Branch("Track_k_m_cov_vy_vz", &track_k_m_cov_vy_vz);

	OutputTree->Branch("Track_k_m_beta", &track_k_m_beta);
	OutputTree->Branch("Track_k_m_beta_err", &track_k_m_beta_err);

	OutputTree->Branch("Track_k_m_filterchi", &track_k_m_filter_chi);
	OutputTree->Branch("Track_k_m_smoothchi", &track_k_m_smooth_chi);
	OutputTree->Branch("Track_k_m_smooth_chi_sum", &track_k_m_smooth_chi_sum);

	OutputTree->Branch("x_estimates_m", &x_estimates_m);
	OutputTree->Branch("y_estimates_m", &y_estimates_m);
	OutputTree->Branch("z_estimates_m", &z_estimates_m);

	OutputTree->Branch("NumTracks_k_m", &NumTracks_k_m);
	OutputTree->Branch("Track_k_m_numHits", &track_k_m_numHits);
	OutputTree->Branch("Track_k_m_hitIndices", &track_hit_k_m_indices);
	OutputTree->Branch("king_move_inds", &king_move_inds);

	OutputTree->Branch("Track_k_m_expected_hit_layer", &track_k_m_expected_hit_layer);
	OutputTree->Branch("Track_k_m_opening_angle", &track_k_m_openingangle);
	OutputTree->Branch("Track_k_m_pdgs", &track_pdgs);
	OutputTree->Branch("Track_k_m_ids", &track_ids);

	OutputTree->Branch("Track_k_m_x_std_scat_per_m", &x_scat);
	OutputTree->Branch("Track_k_m_z_std_scat_per_m", &z_scat);

	//Copy of branches for set beta value
		
	OutputTree->Branch("Vertex_c_b_t", &vertex_c_b_t);
	OutputTree->Branch("Vertex_c_b_x", &vertex_c_b_x);
	OutputTree->Branch("Vertex_c_b_y", &vertex_c_b_y);
	OutputTree->Branch("Vertex_c_b_z", &vertex_c_b_z);
      	OutputTree->Branch("Vertex_c_b_ErrorT", &vertex_t_c_b_error);
      	OutputTree->Branch("Vertex_c_b_ErrorX", &vertex_x_c_b_error);
      	OutputTree->Branch("Vertex_c_b_ErrorY", &vertex_y_c_b_error);
      	OutputTree->Branch("Vertex_c_b_ErrorZ", &vertex_z_c_b_error);
	OutputTree->Branch("vertex_c_b_chi2", &vertex_c_b_chi2_per_dof);
	OutputTree->Branch("Vertex_c_b_trackIndices", &vertex_track_c_b_indices);
	OutputTree->Branch("NumVertices_c_b", &numvertices_c_b, "NumVertices/D");

	OutputTree->Branch("Track_c_b_velX", &track_c_b_vx);
	OutputTree->Branch("Track_c_b_velY", &track_c_b_vy);
	OutputTree->Branch("Track_c_b_velZ", &track_c_b_vz);
	OutputTree->Branch("Track_c_b_x0", &track_c_b_x);
	OutputTree->Branch("Track_c_b_y0", &track_c_b_y);
	OutputTree->Branch("Track_c_b_z0", &track_c_b_z);
	OutputTree->Branch("Track_c_b_t0", &track_c_b_t);

      	OutputTree->Branch("Track_c_b_ErrorT0", &track_c_b_t_error);
      	OutputTree->Branch("Track_c_b_ErrorX0", &track_c_b_x_error);
      	OutputTree->Branch("Track_c_b_ErrorY0", &track_c_b_y_error);
      	OutputTree->Branch("Track_c_b_ErrorZ0", &track_c_b_z_error);
      	OutputTree->Branch("Track_c_b_ErrorVx", &track_c_b_vx_error);
      	OutputTree->Branch("Track_c_b_ErrorVy", &track_c_b_vy_error);
      	OutputTree->Branch("Track_c_b_ErrorVz", &track_c_b_vz_error);

	OutputTree->Branch("Track_c_b_cov_x_t", &track_c_b_cov_x_t);
	OutputTree->Branch("Track_c_b_cov_x_z", &track_c_b_cov_x_z);
	OutputTree->Branch("Track_c_b_cov_x_vx", &track_c_b_cov_x_vx);
	OutputTree->Branch("Track_c_b_cov_x_vy", &track_c_b_cov_x_vy);
	OutputTree->Branch("Track_c_b_cov_x_vz", &track_c_b_cov_x_vz);

	OutputTree->Branch("Track_c_b_cov_t_z", &track_c_b_cov_t_z);
	OutputTree->Branch("Track_c_b_cov_t_vx", &track_c_b_cov_t_vx);
	OutputTree->Branch("Track_c_b_cov_t_vy", &track_c_b_cov_t_vy);
	OutputTree->Branch("Track_c_b_cov_t_vz", &track_c_b_cov_t_vz);

	OutputTree->Branch("Track_c_b_cov_z_vx", &track_c_b_cov_z_vx);
	OutputTree->Branch("Track_c_b_cov_z_vy", &track_c_b_cov_z_vy);
	OutputTree->Branch("Track_c_b_cov_z_vz", &track_c_b_cov_z_vz);

	OutputTree->Branch("Track_c_b_cov_vx_vy", &track_c_b_cov_vx_vy);
	OutputTree->Branch("Track_c_b_cov_vx_vz", &track_c_b_cov_vx_vz);

	OutputTree->Branch("Track_c_b_cov_vy_vz", &track_c_b_cov_vy_vz);

	OutputTree->Branch("Track_c_b_beta", &track_c_b_beta);
	OutputTree->Branch("Track_c_b_beta_err", &track_c_b_beta_err);

	OutputTree->Branch("Track_c_b_filterchi", &track_c_b_filter_chi);
	OutputTree->Branch("Track_c_b_smoothchi", &track_c_b_smooth_chi);
	OutputTree->Branch("Track_c_b_smooth_chi_sum", &track_c_b_smooth_chi_sum);

	OutputTree->Branch("x_estimates_c_b", &x_estimates_c_b);
	OutputTree->Branch("y_estimates_c_b", &y_estimates_c_b);
	OutputTree->Branch("z_estimates_c_b", &z_estimates_c_b);

	OutputTree->Branch("NumTracks_c_b", &NumTracks_c_b);
	OutputTree->Branch("Track_c_b_numHits", &track_c_b_numHits);
	OutputTree->Branch("Track_c_b_hitIndices", &track_hit_c_b_indices);
	OutputTree->Branch("king_move_inds", &c_b_king_move_inds);

	OutputTree->Branch("Track_c_b_expected_hit_layer", &track_c_b_expected_hit_layer);
	OutputTree->Branch("Track_c_b_opening_angle", &track_c_b_openingangle);
	OutputTree->Branch("Track_c_b_pdgs", &track_c_b_pdgs);
	OutputTree->Branch("Track_c_b_ids", &track_c_b_ids);

	OutputTree->Branch("Track_c_b_x_std_scat_per_m", &x_c_b_scat);
	OutputTree->Branch("Track_c_b_z_std_scat_per_m", &z_c_b_scat);

		
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
	std::vector<double> vertex_k_m_chi2_per_dof;
	std::vector<double>	vertex_track_k_m_indices;
	Double_t numvertices_k_m;
	
	//Copy of branches for set beta value
	std::vector<double> q_s_x_m_c_b;
	std::vector<double> q_s_y_m_c_b;
	std::vector<double> q_s_z_m_c_b;
	std::vector<double> vertex_k_f_beta_c_b;
	std::vector<double> vertex_k_s_beta_c_b;
	std::vector<double> vertex_k_s_beta_err_c_b;
	
	std::vector<double>	vertex_c_b_t;
	std::vector<double>	vertex_c_b_x;
	std::vector<double>	vertex_c_b_y;
	std::vector<double>	vertex_c_b_z;
  	std::vector<double> vertex_t_c_b_error;
  	std::vector<double> vertex_x_c_b_error;
  	std::vector<double> vertex_y_c_b_error;
  	std::vector<double> vertex_z_c_b_error;
	std::vector<double> vertex_c_b_chi2_per_dof;
	std::vector<double>	vertex_track_c_b_indices;
	Double_t numvertices_c_b;


  //__Make Track Branches_________________________________________________________________________

	std::vector<double> x_scat;
	std::vector<double> z_scat;

	std::vector<int> track_pdgs;
	std::vector<int> track_ids;
	std::vector<double> track_k_m_openingangle;

	std::vector<double> track_k_m_beta;
	std::vector<double> track_k_m_beta_err;

	std::vector<double> track_k_m_vx;
	std::vector<double> track_k_m_vy;
	std::vector<double> track_k_m_vz;
	std::vector<double> track_k_m_t;
	std::vector<double> track_k_m_x;
	std::vector<double> track_k_m_y;
	std::vector<double> track_k_m_z;

  	std::vector<double> track_k_m_t_error;
  	std::vector<double> track_k_m_x_error;
  	std::vector<double> track_k_m_y_error;
  	std::vector<double> track_k_m_z_error;
  	std::vector<double> track_k_m_vx_error;
  	std::vector<double> track_k_m_vy_error;
  	std::vector<double> track_k_m_vz_error;

	std::vector<double> track_k_m_cov_x_t;
	std::vector<double> track_k_m_cov_x_z;
	std::vector<double> track_k_m_cov_x_vx;
	std::vector<double> track_k_m_cov_x_vy;
	std::vector<double> track_k_m_cov_x_vz;

	std::vector<double> track_k_m_cov_t_z;
	std::vector<double> track_k_m_cov_t_vx;
	std::vector<double> track_k_m_cov_t_vy;
	std::vector<double> track_k_m_cov_t_vz;

	std::vector<double> track_k_m_cov_z_vx;
	std::vector<double> track_k_m_cov_z_vy;
	std::vector<double> track_k_m_cov_z_vz;

	std::vector<double> track_k_m_cov_vx_vy;
	std::vector<double> track_k_m_cov_vx_vz;

	std::vector<double> track_k_m_cov_vy_vz;

	std::vector<double> track_k_m_filter_chi;
	std::vector<double> track_k_m_smooth_chi;
	std::vector<double> track_k_m_smooth_chi_sum;

	int NumTracks_k_m;
	std::vector<double> track_k_m_numHits;
	std::vector<int> track_hit_k_m_indices;

	std::vector<double> x_estimates_m;
	std::vector<double> y_estimates_m;
	std::vector<double> z_estimates_m;

	std::vector<double> track_k_m_expected_hit_layer;
	std::vector<int> king_move_inds;
	
	//Copy of branches for set beta value
	
	std::vector<double> x_c_b_scat;
	std::vector<double> z_c_b_scat;

	std::vector<int> track_c_b_pdgs;
	std::vector<int> track_c_b_ids;
	std::vector<double> track_c_b_openingangle;

	std::vector<double> track_c_b_beta;
	std::vector<double> track_c_b_beta_err;

	std::vector<double> track_c_b_vx;
	std::vector<double> track_c_b_vy;
	std::vector<double> track_c_b_vz;
	std::vector<double> track_c_b_t;
	std::vector<double> track_c_b_x;
	std::vector<double> track_c_b_y;
	std::vector<double> track_c_b_z;

  	std::vector<double> track_c_b_t_error;
  	std::vector<double> track_c_b_x_error;
  	std::vector<double> track_c_b_y_error;
  	std::vector<double> track_c_b_z_error;
  	std::vector<double> track_c_b_vx_error;
  	std::vector<double> track_c_b_vy_error;
  	std::vector<double> track_c_b_vz_error;

	std::vector<double> track_c_b_cov_x_t;
	std::vector<double> track_c_b_cov_x_z;
	std::vector<double> track_c_b_cov_x_vx;
	std::vector<double> track_c_b_cov_x_vy;
	std::vector<double> track_c_b_cov_x_vz;

	std::vector<double> track_c_b_cov_t_z;
	std::vector<double> track_c_b_cov_t_vx;
	std::vector<double> track_c_b_cov_t_vy;
	std::vector<double> track_c_b_cov_t_vz;

	std::vector<double> track_c_b_cov_z_vx;
	std::vector<double> track_c_b_cov_z_vy;
	std::vector<double> track_c_b_cov_z_vz;

	std::vector<double> track_c_b_cov_vx_vy;
	std::vector<double> track_c_b_cov_vx_vz;

	std::vector<double> track_c_b_cov_vy_vz;

	std::vector<double> track_c_b_filter_chi;
	std::vector<double> track_c_b_smooth_chi;
	std::vector<double> track_c_b_smooth_chi_sum;

	int NumTracks_c_b;
	std::vector<double> track_c_b_numHits;
	std::vector<int> track_hit_c_b_indices;

	std::vector<double> x_estimates_c_b;
	std::vector<double> y_estimates_c_b;
	std::vector<double> z_estimates_c_b;

	std::vector<double> track_c_b_expected_hit_layer;
	std::vector<int> c_b_king_move_inds;


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
    long long int digi_seed;

}; //class TreeHandler

template<class digi_hit>
void TreeHandler::ExportDigis(std::vector<digi_hit*> digi_list, long long int seed){
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
      
      digi_seed = seed;

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
void TreeHandler::ExportTracks_k_m(std::vector<Track*> track_list){ 
  // Data pushed to Tree for all the Tracks (that survived the merging algorithm)

	track_k_m_vx.clear();
	track_k_m_vy.clear();
	track_k_m_vz.clear();
	track_k_m_x.clear();
	track_k_m_y.clear();
	track_k_m_z.clear();
	track_k_m_t.clear();

	track_k_m_vx_error.clear();
	track_k_m_vy_error.clear();
	track_k_m_vz_error.clear();
	track_k_m_x_error.clear();
	track_k_m_y_error.clear();
	track_k_m_z_error.clear();
	track_k_m_t_error.clear();
	
	x_estimates_m.clear();
	y_estimates_m.clear();
	z_estimates_m.clear();

	track_k_m_cov_x_t.clear();
	track_k_m_cov_x_z.clear();
	track_k_m_cov_x_vx.clear();
	track_k_m_cov_x_vy.clear();
	track_k_m_cov_x_vz.clear();

	track_k_m_cov_t_z.clear();
	track_k_m_cov_t_vx.clear();
	track_k_m_cov_t_vy.clear();
	track_k_m_cov_t_vz.clear();

	track_k_m_cov_z_vx.clear();
	track_k_m_cov_z_vy.clear();
	track_k_m_cov_z_vz.clear();

	track_k_m_cov_vx_vy.clear();
	track_k_m_cov_vx_vz.clear();

	track_k_m_cov_vy_vz.clear();

  track_hit_k_m_indices.clear();
  track_k_m_expected_hit_layer.clear();
	king_move_inds.clear();
 
  track_k_m_filter_chi.clear();
  track_k_m_smooth_chi.clear();
	track_k_m_smooth_chi_sum.clear();
 
  NumTracks_k_m = track_list.size();
  track_k_m_numHits.clear();

  track_k_m_openingangle.clear();
  track_ids.clear();
  track_pdgs.clear();
  
  x_scat.clear();
  z_scat.clear();
  
  track_k_m_beta.clear();
  track_k_m_beta_err.clear();
  
  // Opening angles among all pairs of tracks
  for (int i = 0; i < NumTracks_k_m; i++) {
  		for (int j = i+1; j < NumTracks_k_m; j++) {
  		auto tr1 =  track_list[i];
  		auto tr2 =  track_list[j];
  		track_k_m_openingangle.push_back(tr1->direction()^tr2->direction());
  	}
  }

  for (auto tr : track_list){
      // Number of Tracks Reconstructed
      track_k_m_numHits.push_back( (tr->hits).size() );

     /*
      if ((tr->hits).size()  > 8 ) {
	std::cout << "num_hits " << (tr->hits).size() << std::endl;
\	for (auto hit : tr->hits) {
		std::cout << " hit y: " << hit->y;
	}
	std::cout << std::endl;
	for (auto hit : tr->hits) {
		std::cout << " hit z: " << hit->z;
	}
	std::cout << std::endl;
      }
      */

      // Push velocity and position of lowest (in y) Kalman estimate for each track
    	track_k_m_vx.push_back(tr->vx);
    	track_k_m_vy.push_back(tr->vy);
    	track_k_m_vz.push_back(tr->vz);
     
    	track_k_m_t.push_back(tr->t0);
    	track_k_m_x.push_back(tr->x0);
    	track_k_m_y.push_back(tr->y0);
    	track_k_m_z.push_back(tr->z0);

    	track_k_m_vx_error.push_back(tr->evx);
    	track_k_m_vy_error.push_back(tr->evy);
    	track_k_m_vz_error.push_back(tr->evz);
     
    	track_k_m_t_error.push_back(tr->et0);
    	track_k_m_x_error.push_back(tr->ex0);
    	track_k_m_y_error.push_back(tr->ey0);
    	track_k_m_z_error.push_back(tr->ez0);
    
	track_k_m_cov_x_t.push_back(tr->P_s(0,1));
	track_k_m_cov_x_z.push_back(tr->P_s(0,2));
	track_k_m_cov_x_vx.push_back(tr->P_s(0,3));
	track_k_m_cov_x_vy.push_back(tr->P_s(0,4));
	track_k_m_cov_x_vz.push_back(tr->P_s(0,5));

	track_k_m_cov_t_z.push_back(tr->P_s(1,2));
	track_k_m_cov_t_vx.push_back(tr->P_s(1,3));
	track_k_m_cov_t_vy.push_back(tr->P_s(1,4));
	track_k_m_cov_t_vz.push_back(tr->P_s(1,5));

	track_k_m_cov_z_vx.push_back(tr->P_s(2,3));
	track_k_m_cov_z_vy.push_back(tr->P_s(2,4));
	track_k_m_cov_z_vz.push_back(tr->P_s(2,5));

	track_k_m_cov_vx_vy.push_back(tr->P_s(3,4));
	track_k_m_cov_vx_vz.push_back(tr->P_s(3,5));

	track_k_m_cov_vy_vz.push_back(tr->P_s(4,5));




      // Push various Chis to Tree
      for (auto chi : tr->chi_f) track_k_m_filter_chi.push_back( static_cast<double>(chi) );
      for (auto chi : tr->chi_s) track_k_m_smooth_chi.push_back( static_cast<double>(chi) );
      track_k_m_filter_chi.push_back(-1); // Raw chi increments at each kalman step of the final fit
      track_k_m_smooth_chi.push_back(-1); 

      double chi_sum = 0;
			for (auto chi : tr->chi_s) {
				chi_sum += chi;
			}
			chi_sum = chi_sum / (4.0*tr->chi_s.size() - 6.0);
			track_k_m_smooth_chi_sum.push_back(chi_sum); // final chi per ndof for the track

      // Indices of various sorts
    	for (auto ind : tr->king_move_inds) {
    		king_move_inds.push_back(ind);
    	}
    	if (tr->king_move_inds.size() != 0) king_move_inds.push_back(-1); // indices of hits removed by king moves algorithm
                                                                         // see presentation slides in Doc for description

    	for (auto hit : tr->hits) {
    		track_hit_k_m_indices.push_back(hit->index);
    	}
    	track_hit_k_m_indices.push_back(-1.); // indices of digi hits included in the tracks

      // Indices of layers where tracks are expected to have travelled through
    	for (auto exp_layer : tr->expected_layers) {
    		track_k_m_expected_hit_layer.push_back( static_cast<double>(exp_layer) );
    	}
    	track_k_m_expected_hit_layer.push_back(-1);

      // Scattering per m calculated for covariance parallel to detector planes
      for (auto std : tr->x_scats) x_scat.push_back(std);
			for (auto std : tr->z_scats) z_scat.push_back(std);

      // List of all kalman best estimates of hit positions for each track
    	for (auto est : tr->estimate_list) {
    		x_estimates_m.push_back(est[0]);
    		y_estimates_m.push_back(est[1]);
    		z_estimates_m.push_back(est[2]);
    	}
    	x_estimates_m.push_back(-1.0);
    	y_estimates_m.push_back(-1.0);
    	z_estimates_m.push_back(-1.0);

      // Track and Pdg Ids of hits
      for (auto hit : tr->hits) {
				track_ids.push_back(hit->min_track_id);
				track_pdgs.push_back(hit->pdg);
			}
			track_ids.push_back(-1);
			track_pdgs.push_back(-1);

      // Track Betas and Error in Beta
      double v = std::sqrt(std::pow(tr->vx,2) + std::pow(tr->vy,2) + std::pow(tr->vz,2));
			track_k_m_beta.push_back(v / constants::c);

			Eigen::MatrixXd R = tr->P_s;
			Eigen::MatrixXd D(3,3);
			D << R(3,3), R(3,4), R(3,5),
		             R(4,3), R(4,4), R(4,5),
       			     R(5,3), R(5,4), R(5,5);
			Eigen::VectorXd q(3);
			q << tr->vx, tr->vy, tr->vz;
			track_k_m_beta_err.push_back((q.transpose() * D * q)(0) / (v * v));

      }

}

template<class Track>
void TreeHandler::ExportTracks_c_b(std::vector<Track*> track_list){ 
  // Data pushed to Tree for all the Tracks (that survived the merging algorithm)

	track_c_b_vx.clear();
	track_c_b_vy.clear();
	track_c_b_vz.clear();
	track_c_b_x.clear();
	track_c_b_y.clear();
	track_c_b_z.clear();
	track_c_b_t.clear();

	track_c_b_vx_error.clear();
	track_c_b_vy_error.clear();
	track_c_b_vz_error.clear();
	track_c_b_x_error.clear();
	track_c_b_y_error.clear();
	track_c_b_z_error.clear();
	track_c_b_t_error.clear();
	
	x_estimates_c_b.clear();
	y_estimates_c_b.clear();
	z_estimates_c_b.clear();

	track_c_b_cov_x_t.clear();
	track_c_b_cov_x_z.clear();
	track_c_b_cov_x_vx.clear();
	track_c_b_cov_x_vy.clear();
	track_c_b_cov_x_vz.clear();

	track_c_b_cov_t_z.clear();
	track_c_b_cov_t_vx.clear();
	track_c_b_cov_t_vy.clear();
	track_c_b_cov_t_vz.clear();

	track_c_b_cov_z_vx.clear();
	track_c_b_cov_z_vy.clear();
	track_c_b_cov_z_vz.clear();

	track_c_b_cov_vx_vy.clear();
	track_c_b_cov_vx_vz.clear();

	track_c_b_cov_vy_vz.clear();

  track_hit_c_b_indices.clear();
  track_c_b_expected_hit_layer.clear();
	c_b_king_move_inds.clear();
 
  track_c_b_filter_chi.clear();
  track_c_b_smooth_chi.clear();
	track_c_b_smooth_chi_sum.clear();
 
  NumTracks_c_b = track_list.size();
  track_c_b_numHits.clear();

  track_c_b_openingangle.clear();
  track_c_b_ids.clear();
  track_c_b_pdgs.clear();
  
  x_c_b_scat.clear();
  z_c_b_scat.clear();
  
  track_c_b_beta.clear();
  track_c_b_beta_err.clear();
  
  // Opening angles among all pairs of tracks
  for (int i = 0; i < NumTracks_c_b; i++) {
  		for (int j = i+1; j < NumTracks_c_b; j++) {
  		auto tr1 =  track_list[i];
  		auto tr2 =  track_list[j];
  		track_c_b_openingangle.push_back(tr1->direction()^tr2->direction());
  	}
  }

  for (auto tr : track_list){
      // Number of Tracks Reconstructed
      track_c_b_numHits.push_back( (tr->hits).size() );

     /*
      if ((tr->hits).size()  > 8 ) {
	std::cout << "c_b num_hits " << (tr->hits).size() << std::endl;
\	for (auto hit : tr->hits) {
		std::cout << "c_b  hit y: " << hit->y;
	}
	std::cout << std::endl;
	for (auto hit : tr->hits) {
		std::cout << "c_b hit z: " << hit->z;
	}
	std::cout << std::endl;
      }
      */

      // Push velocity and position of lowest (in y) Kalman estimate for each track
    	track_c_b_vx.push_back(tr->vx);
    	track_c_b_vy.push_back(tr->vy);
    	track_c_b_vz.push_back(tr->vz);
     
    	track_c_b_t.push_back(tr->t0);
    	track_c_b_x.push_back(tr->x0);
    	track_c_b_y.push_back(tr->y0);
    	track_c_b_z.push_back(tr->z0);

    	track_c_b_vx_error.push_back(tr->evx);
    	track_c_b_vy_error.push_back(tr->evy);
    	track_c_b_vz_error.push_back(tr->evz);
     
    	track_c_b_t_error.push_back(tr->et0);
    	track_c_b_x_error.push_back(tr->ex0);
    	track_c_b_y_error.push_back(tr->ey0);
    	track_c_b_z_error.push_back(tr->ez0);
    
	track_c_b_cov_x_t.push_back(tr->P_s(0,1));
	track_c_b_cov_x_z.push_back(tr->P_s(0,2));
	track_c_b_cov_x_vx.push_back(tr->P_s(0,3));
	track_c_b_cov_x_vy.push_back(tr->P_s(0,4));
	track_c_b_cov_x_vz.push_back(tr->P_s(0,5));

	track_c_b_cov_t_z.push_back(tr->P_s(1,2));
	track_c_b_cov_t_vx.push_back(tr->P_s(1,3));
	track_c_b_cov_t_vy.push_back(tr->P_s(1,4));
	track_c_b_cov_t_vz.push_back(tr->P_s(1,5));

	track_c_b_cov_z_vx.push_back(tr->P_s(2,3));
	track_c_b_cov_z_vy.push_back(tr->P_s(2,4));
	track_c_b_cov_z_vz.push_back(tr->P_s(2,5));

	track_c_b_cov_vx_vy.push_back(tr->P_s(3,4));
	track_c_b_cov_vx_vz.push_back(tr->P_s(3,5));

	track_c_b_cov_vy_vz.push_back(tr->P_s(4,5));




      // Push various Chis to Tree
      for (auto chi : tr->chi_f) track_c_b_filter_chi.push_back( static_cast<double>(chi) );
      for (auto chi : tr->chi_s) track_c_b_smooth_chi.push_back( static_cast<double>(chi) );
      track_c_b_filter_chi.push_back(-1); // Raw chi increments at each kalman step of the final fit
      track_c_b_smooth_chi.push_back(-1); 

      double chi_sum = 0;
			for (auto chi : tr->chi_s) {
				chi_sum += chi;
			}
			chi_sum = chi_sum / (4.0*tr->chi_s.size() - 6.0);
			track_c_b_smooth_chi_sum.push_back(chi_sum); // final chi per ndof for the track

      // Indices of various sorts
	//TODO: make sure the king moves are accurately done here with c_b added in
    	for (auto ind : tr->king_move_inds) {
    		c_b_king_move_inds.push_back(ind);
    	}
    	if (tr->c_b_king_move_inds.size() != 0) c_b_king_move_inds.push_back(-1); // indices of hits removed by king moves algorithm
                                                                         // see presentation slides in Doc for description

    	for (auto hit : tr->hits) {
    		track_hit_c_b_indices.push_back(hit->index);
    	}
    	track_hit_c_b_indices.push_back(-1.); // indices of digi hits included in the tracks

      // Indices of layers where tracks are expected to have travelled through
    	for (auto exp_layer : tr->expected_layers) {
    		track_c_b_expected_hit_layer.push_back( static_cast<double>(exp_layer) );
    	}
    	track_c_b_expected_hit_layer.push_back(-1);

      // Scattering per m calculated for covariance parallel to detector planes
      for (auto std : tr->x_scats) x_c_b_scat.push_back(std);
			for (auto std : tr->z_scats) z_c_b_scat.push_back(std);

      // List of all kalman best estimates of hit positions for each track
    	for (auto est : tr->estimate_list) {
    		x_estimates_c_b.push_back(est[0]);
    		y_estimates_c_b.push_back(est[1]);
    		z_estimates_c_b.push_back(est[2]);
    	}
    	x_estimates_c_b.push_back(-1.0);
    	y_estimates_c_b.push_back(-1.0);
    	z_estimates_c_b.push_back(-1.0);

      // Track and Pdg Ids of hits
      for (auto hit : tr->hits) {
				track_ids.push_back(hit->min_track_id);
				track_pdgs.push_back(hit->pdg);
			}
			track_c_b_ids.push_back(-1);
			track_c_b_pdgs.push_back(-1);

      // Track Betas and Error in Beta
      double v = std::sqrt(std::pow(tr->vx,2) + std::pow(tr->vy,2) + std::pow(tr->vz,2));
			track_c_b_beta.push_back(v / constants::c);

			Eigen::MatrixXd R = tr->P_s;
			Eigen::MatrixXd D(3,3);
			D << R(3,3), R(3,4), R(3,5),
		             R(4,3), R(4,4), R(4,5),
       			     R(5,3), R(5,4), R(5,5);
			Eigen::VectorXd q(3);
			q << tr->vx, tr->vy, tr->vz;
			track_c_b_beta_err.push_back((q.transpose() * D * q)(0) / (v * v));

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
	vertex_k_m_chi2_per_dof.clear();
	numvertices_k_m = vertices.size();

	for (auto v : vertices) {

		vertex_k_m_t.push_back(v->t);
		vertex_k_m_x.push_back(v->x);
		vertex_k_m_y.push_back(v->y);
		vertex_k_m_z.push_back(v->z);

		vertex_x_k_m_error.push_back(sqrt(v->CovMatrix()[0][0]));
		vertex_y_k_m_error.push_back(sqrt(v->CovMatrix()[1][1]));
		vertex_z_k_m_error.push_back(sqrt(v->CovMatrix()[2][2]));
		vertex_t_k_m_error.push_back(sqrt(v->CovMatrix()[3][3]));

		vertex_k_m_chi2_per_dof.push_back(v->merit());

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

	        vertex_track_k_m_indices.push_back(-1);
	}

}

template<typename vertex>
void TreeHandler::ExportVertices_c_b(std::vector<vertex*> vertices){

	vertex_c_b_t.clear();
	vertex_c_b_x.clear();
	vertex_c_b_y.clear();
	vertex_c_b_z.clear();
	vertex_track_c_b_indices.clear();

	q_s_x_m_c_b.clear();
	q_s_y_m_c_b.clear();
	q_s_z_m_c_b.clear();
	vertex_k_f_beta_c_b.clear();
	vertex_k_s_beta_c_b.clear();
	vertex_k_s_beta_err_c_b.clear();

	vertex_t_c_b_error.clear();
	vertex_x_c_b_error.clear();
	vertex_y_c_b_error.clear();
	vertex_z_c_b_error.clear();
	vertex_c_b_chi2_per_dof.clear();
	numvertices_c_b = vertices.size();

	for (auto v : vertices) {

		vertex_c_b_t.push_back(v->t);
		vertex_c_b_x.push_back(v->x);
		vertex_c_b_y.push_back(v->y);
		vertex_c_b_z.push_back(v->z);

		vertex_x_c_b_error.push_back(sqrt(v->CovMatrix()[0][0]));
		vertex_y_c_b_error.push_back(sqrt(v->CovMatrix()[1][1]));
		vertex_z_c_b_error.push_back(sqrt(v->CovMatrix()[2][2]));
		vertex_t_c_b_error.push_back(sqrt(v->CovMatrix()[3][3]));

		vertex_c_b_chi2_per_dof.push_back(v->merit());

		for (auto q_f : v->q_f) vertex_k_f_beta_c_b.push_back(q_f.norm() / constants::c);
		vertex_k_f_beta_c_b.push_back(-1.0);

		for (int i=0; i < v->q_s.size(); i++) {
			q_s_x_m_c_b.push_back(v->q_s[i][0]);
			q_s_y_m_c_b.push_back(v->q_s[i][1]);
			q_s_z_m_c_b.push_back(v->q_s[i][2]);

			//double pull = (v->q_s[i].norm() - constants::c) * v->q_s[i].squaredNorm() /
			//	(v->q_s[i].transpose() * v->D_s[i] * v->q_s[i])(0);
			//vertex_k_beta_c_b.push_back(pull);

			double vertex_k_beta_err = (v->q_s[i].transpose() * v->D_s[i] * v->q_s[i])(0) / v->q_s[i].squaredNorm();

			vertex_k_s_beta_err_c_b.push_back(vertex_k_beta_err);

			vertex_k_s_beta_c_b.push_back((v->q_s[i]).norm() / constants::c);
		}
		q_s_x_m_c_b.push_back(-1.0);
		q_s_y_m_c_b.push_back(-1.0);
		q_s_z_m_c_b.push_back(-1.0);
		vertex_k_s_beta_c_b.push_back(-1.0);

		for (auto tr_index : v->track_indices){
        		vertex_track_c_b_indices.push_back(tr_index);
        	}

	        vertex_track_c_b_indices.push_back(-1);
	}

}



#endif
