#include <iostream>
#include "RunManager.hh"
#include "TString.h"
#include "util.hh"


int main(int argc, char *argv[]){

	if (argc != 3) {
		std::cout << "Need 2 Arguments! \n First argument: input_file_name (or dir) \n Second argument: output_dir_name" << std::endl;
		return 0;
	}

	std::vector<TString> files;

	if ( io::is_directory( argv[1] ) ) {
		files = io::ProcessDirectory(argv[1], "");
	} else {
		files = {TString(argv[1])};
	}
	
	TString outdir = TString(argv[2]);

	RunManager RM;
	RM.SetOutputFile(outdir);

	for (auto f : files ){
		std::cout << f << std::endl;
		RM.SetInputFile(f);
		RM.StartTracking();
	} 

	

		


	return 0;

} //main