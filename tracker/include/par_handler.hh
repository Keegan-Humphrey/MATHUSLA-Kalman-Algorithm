#include <map>

#pragma once

class ParHandler {

public:
	std::map<std::string, double> par_map;

	void Handle();

	bool file_opened;
private:
	void GetPars();

};


