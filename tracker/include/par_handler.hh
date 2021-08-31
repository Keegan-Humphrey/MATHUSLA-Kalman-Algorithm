#include <map>

#pragma once

class ParHandler {

public:
	std::map<std::string, double> par_map;

	void Handle();

private:
	void GetPars();

};


