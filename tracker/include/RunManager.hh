#include <TString.h>
#include "TrackFinder.hh"
#include "TreeHandler.hh"
#include "Digitizer.hh"
#include "globals.hh"
#include "VertexFinder.hh"


#ifndef RUN_MANAGER_DEFINE
#define RUN_MANAGER_DEFINE

class RunManager{


public:

	// called in main function to start algorithm
	int StartTracking();
	void SetInputFile(TString name){_InputFile_Name = name; fileNumber++;}
	void SetOutputFile(TString name){_OutputFile_Name = name;}
	TString OutFileName(){ 
		std::ostringstream strs;
		strs << fileNumber;
		return _OutputFile_Name + TString("/stat") + TString(strs.str()) + TString(".root"); 
	}

	unsigned long long TotalEventsProcessed = 0;

	RunManager(){
		
		_digitizer = new Digitizer();
		_tracker = new TrackFinder();
		_vertexer = new VertexFinder();
	} 

	~RunManager(){
		delete _digitizer;
		delete _tracker;
		delete _vertexer;
	}


private:

	//DATA AND TRACKING VARIABLES
	int fileNumber = -1;
	int LoadEvent(int);
	TreeHandler* TH;

	Digitizer* _digitizer;
	TrackFinder* _tracker;
	VertexFinder* _vertexer;

	//DATA IO NAMES AND FILES
	TString _InputFile_Name;
	TString _OutputFile_Name;

	TString _InputTree_Name = TString("box_run");
	TString _OutputTree_Name = TString("integral_tree");



};


#endif