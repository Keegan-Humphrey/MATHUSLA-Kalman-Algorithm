#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include "par_handler.hh"
//#include <filesystem>



void ParHandler::Handle() {

        GetPars();

}


void ParHandler::GetPars() {

//	std::filesystem::path par_path;

	std::string par_path = std::string(__FILE__);

	std::string par_card_str = par_path.substr(0, par_path.size()-14)+"/../run/par_card.txt";
	//par_card_str.st = myString.substr(0, myString.size()-1);

	std::cout << par_card_str << std::endl;
        std::ifstream infile(par_card_str);

//        std::ifstream infile("../run/par_card.txt");
//        std::ifstream infile("/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/tracker/run/par_card.txt");

	if (infile.is_open()) {
		file_opened = true;

	        std::string line;
        	while (std::getline(infile, line))
        	{
                	std::istringstream iss(line);
                	std::string par;
                	double value;

	                if (!(iss >> par >> value)) {
	//			std::cout << "Sorry I couldn't find the parameters I need to run." << std::endl;
				continue; } // error

			//std::cout << par << " " << value << std::endl;

	                par_map[par] = value;
        	}
	}
	else {
		file_opened = false;

		std::cout << "Error opening parameter card file" << std::endl;
	}
}

