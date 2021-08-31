#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include "par_handler.hh"




void ParHandler::Handle() {

        GetPars();

}


void ParHandler::GetPars() {

        std::ifstream infile("../run/par_card.txt");

        std::string line;
        while (std::getline(infile, line))
        {
                std::istringstream iss(line);
                std::string par;
                double value;

                if (!(iss >> par >> value)) { continue; } // error

                par_map[par] = value;
        }

}

