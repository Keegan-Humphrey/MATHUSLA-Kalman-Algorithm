#include "globals.hh"
#include "LinearAlgebra.hh"
#include <cmath>
#include <iostream>

#ifndef GEOMETRY_HH
#define GEOMETRY_HH

class Layer  {
public:
	int index;
	double min;
	double max;
	double center;
	int long_direction_index;
	int short_direction_index;

	//the uncertainty along each directiion in the layer
	std::vector<double> widths(){
		if ( (index % 2 ) == 0){
			return { detector::scintillator_width, detector::scintillator_length };
		}

		return {  detector::scintillator_length, detector::scintillator_width  };
	}
	std::vector<double> uncertainty(){
		if ( (index % 2 ) == 0){
			return { detector::scintillator_width/sqrt(12.), detector::scintillator_thickness/sqrt(12.),
                detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2) };
		}

		return {  detector::time_resolution * (constants::c/constants::optic_fiber_n)/ sqrt(2), (max-min)/sqrt(12.), detector::scintillator_width/sqrt(12.) };

	}

	Layer(int _index, double _min, double _max){
		index = _index;
		min = _min;
		max = _max;
		center = (min + max)/2.0;

		if (index % 2 == 0){
			short_direction_index = 0;
			long_direction_index = 1;
		} else {
			short_direction_index = 1;
			long_direction_index = 0;
		}
	}

	bool in_layer(double y){
		if (y > min and y < max) return true;
		return false;
	 }

};

class Module{
public:
	double cx, cy, cz; //center point
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	int index;
	std::vector<Layer*> layers;
	void SetLayers(std::vector<Layer*> _layers){layers = _layers;}
	void SetIndex(int _index){index = _index;}

	std::vector<double> LocalPosition(double x, double y, double z){
		return {x-cx, y-cy, z-cz};
	}


	Module(double _xmin, double _xmax, double _ymin, double _ymax, double _zmin, double _zmax){
		xmin = _xmin;
		xmax = _xmax;
		ymin = _ymin;
		ymax = _ymax;
		zmin = _zmin;
		zmax = _zmax;

		cx = (xmin + xmax)/2.00;
		cy = (ymin + ymax)/2.00;
		cz = (zmin + zmax)/2.00;

	}

	bool is_inside(double x, double y, double z){
		if ( !(x > xmin and x < xmax ) ) return false;
		if ( !(y > ymin and y < ymax ) ) return false;
		if ( !(z > zmin and z < zmax ) ) return false;
		return true;
	}
};

class Floor {
public:
	double cx = (detector::x_min + detector::x_max)/2.0;
	double cz = (detector::z_min + detector::z_max)/2.0;
	double cy_0 = (detector::LAYERS_Y[0][1] + detector::LAYERS_Y[0][0])/2.0;
	double cy_1 = (detector::LAYERS_Y[1][1] + detector::LAYERS_Y[1][0])/2.0;
	double cy_2 = (detector::LAYERS_Y[2][1] + detector::LAYERS_Y[2][0])/2.0;
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	double x_width = detector::floor_x_width;
	double z_width = detector::floor_z_width;


	std::vector<int> GetFloorIndex(double x, double y, double z){
		double local_x = x - cx;
		//double local_y = y - cy;
		double local_z = z - cz;

		return { (int)(std::floor( local_x/x_width )),  (int)(std::floor( local_z/z_width ))   };

	}

	template<typename _detID_>
	std::vector<double> GetCenter(_detID_ id){
		if (!id.isFloorElement) {
			std::cout << "Not floor element!" << std::endl;
			return {};
		}

		double x_local = x_width*(static_cast<double>(id.xIndex) + 0.5);
		double z_local = z_width*(static_cast<double>(id.zIndex) + 0.5);

		if (id.layerIndex == 0) {
			return {cx + x_local, cy_0, cz + z_local};
		} else if (id.layerIndex == 1){
			return {cx + x_local, cy_1, cz + z_local};
		} else if (id.layerIndex == 2){
			return {cx + x_local, cy_2, cz + z_local};
		}

	}

	std::vector<double> uncertainty(int y_index){
		if (y_index == 0) {
			return {x_width/sqrt(12.0), detector::scintillator_thickness/sqrt(12.), z_width/sqrt(12)};
		} else if (y_index == 1) {
			return {x_width/sqrt(12.0), detector::scintillator_thickness/sqrt(12.), z_width/sqrt(12)};
		} else if (y_index == 2) {
			return {x_width/sqrt(12.0), detector::scintillator_thickness/sqrt(12.), z_width/sqrt(12)};
		}
	}
};

class Wall {
public:
	double cx = (detector::x_min + detector::x_max)/2.0;
	double cz = (detector::z_min - detector::wall_gap - detector::scintillator_height/2.0);
	double cy = (detector::y_min + detector::wall_height/2.0);
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	double x_width = detector::wall_x_width;
	double y_width = detector::wall_y_width;


	std::vector<int> GetWallIndex(double x, double y, double z){
		double local_x = x - cx;
		double local_y = y - cy;
		//double local_z = z - cz;

		return { (int)(std::floor( local_x/x_width )),  (int)(std::floor( local_y/y_width ))   };

	}

	template<typename _detID_>
	std::vector<double> GetCenter(_detID_ id){
		if (!id.isWallElement) {
			std::cout << "Not wall element!" << std::endl;
			return {};
		}

		double x_local = x_width*(static_cast<double>(id.xIndex) + 0.5);
		double y_local = y_width*(static_cast<double>(id.wall_yIndex) + 0.5);

		return { cx + x_local, cy + y_local, cz }; // z is getting fixed here for hit < 7000

	}

	std::vector<double> uncertainty(int z_index = 0){
        return { x_width/sqrt(12.0), y_width/sqrt(12.0), detector::scintillator_thickness/sqrt(12.0) };
	}
};



class detID{
public:
	int moduleIndex;
	int layerIndex;
	int xIndex;
	int zIndex;
	bool _null = false;
	bool isFloorElement = false;
    bool isWallElement = false;
    int wall_yIndex = 0;

	detID(){_null = true;}

	void Print(){

		std::cout << "***************Printing DetID*************" << std::endl;
		if (_null){
			std::cout << "NOT IN KNOWN DETECTOR ELEMENT" << std::endl;
			return;
		}
		std::cout << "Module Index: " << moduleIndex << std::endl;
		std::cout << "Layer Index: " << layerIndex << std::endl;
		std::cout << "x Index: " << xIndex << std::endl;
		std::cout << "z Index: " << zIndex << std::endl;
	}

	detID(int _module_index, int _layer_index, int _x_index, int _z_index, bool _isFloorElement = false,
            bool _isWallElement = false, int _wall_y_index = 0){
		moduleIndex = _module_index;
		layerIndex = _layer_index;
		xIndex = _x_index;
		zIndex = _z_index;
		_null = false; // could be ==? used to be, thought it was weird, took it out.
		isFloorElement = _isFloorElement;
        isWallElement = _isWallElement;
        wall_yIndex = _wall_y_index;
	}

	bool IsNull(){return _null;}

	bool operator==(const detID &detID2){
		if (_null) return false;
		if (isFloorElement != detID2.isFloorElement) return false;
        if (isWallElement != detID2.isWallElement) return false;
        
		if (isFloorElement){
			if (xIndex != detID2.xIndex) return false;
			if (zIndex != detID2.zIndex) return false;
			return true;
		}
        if (isWallElement) {
            if (xIndex != detID2.xIndex) return false;
            if (wall_yIndex != detID2.wall_yIndex) return false;
            return true;
        }

		if (moduleIndex != detID2.moduleIndex) return false;
		if (layerIndex != detID2.layerIndex) return false;
		if (xIndex != detID2.xIndex) return false;
		if (zIndex != detID2.zIndex) return false;
		return true;
	}
};

class Geometry{
public:

	std::vector<Module*> module_list{};
	std::vector<Layer*> layer_list{};

	//Extra geometry elements here
	Floor _floor;
    Wall _wall;

	Geometry(){

		for (int _index = 0; _index < detector::n_layers; _index++ ){
			layer_list.push_back(new Layer(_index, detector::LAYERS_Y[_index][0], detector::LAYERS_Y[_index][1]));

		}


		//CONSTRUCT MODULE LIST
		auto min_y = detector::LAYERS_Y[0][0];

		auto max_y = detector::LAYERS_Y[detector::n_layers-1][1];

		for (int _index = 0; _index < detector::n_modules; _index++){
			int _x = _index % 10;
			int _z = (_index - _x)/10;
			module_list.push_back(new Module(detector::MODULE_X[_x][0], detector::MODULE_X[_x][1], min_y, max_y, detector::MODULE_Z[_z][0], detector::MODULE_Z[_z][1]));
			module_list[_index]->SetIndex(_index);
		}

		_floor = Floor();



	} //Geometry Constructor

	~Geometry(){
		for (auto p : module_list){delete p;}
		for (auto p : layer_list){delete p;}
	}


	detID GetDetIDFloor(double x, double y, double z){

		int module_index = -1;
		bool isFloorElement = true;

		if (y < detector::LAYERS_Y[0][1] && y > detector::LAYERS_Y[0][0]) {
			int layer_number = 0;
			std::vector<int> floor_indices = _floor.GetFloorIndex(x, y, z);
			return detID(module_index, layer_number, floor_indices[0], floor_indices[1], isFloorElement);
		} else if (y < detector::LAYERS_Y[1][1] && y > detector::LAYERS_Y[1][0]){
			int layer_number = 1;
			std::vector<int> floor_indices = _floor.GetFloorIndex(x, y, z);
			return detID(module_index, layer_number, floor_indices[0], floor_indices[1], isFloorElement);
		} else if (y < detector::LAYERS_Y[2][1] && y > detector::LAYERS_Y[2][0]){
			int layer_number = 2;
			std::vector<int> floor_indices = _floor.GetFloorIndex(x, y, z);
			return detID(module_index, layer_number, floor_indices[0], floor_indices[1], isFloorElement);
		}

	}

    detID GetDetIDWall(double x, double y, double z) {
        
        int module_index = -1;
        bool isWallElement = true;
        int layer_number = 0;

        std::vector<int> wall_indices = _wall.GetWallIndex(x, y, z);
        return detID(module_index, layer_number, wall_indices[0], 0, false, isWallElement, wall_indices[1]);   
           
    }

	detID GetDetID(double x, double y, double z){
		int module_index = -1;
		int layer_number = -1;
		std::vector<double> layer_widths = {-1, -1};

		//FINDING LAYER

		for (auto layer : layer_list){
			if (layer->in_layer(y)){
				layer_number = layer->index;
				layer_widths = layer->widths();
				break;
			}
		}

//        if (z < 7000) return GetDetIDWall(x, y, z);
        if (z < detector::z_min) return GetDetIDWall(x, y, z);

		if (layer_number == -1){
			return detID();}

		if (layer_number < 3) return GetDetIDFloor(x, y, z);

		//FINDING MODULE

		for (auto module : module_list){

			if (module->is_inside(x, y, z)){
				module_index = module->index;
				break;
			}
		}

		if (module_index == -1){
			return detID();}



		std::vector<double> local_position = (module_list[module_index])->LocalPosition(x, y, z);
		layer_widths = layer_list[layer_number]->widths();

		int x_index = static_cast<size_t>(std::floor(local_position[0]/ layer_widths[0]));
		int z_index = static_cast<size_t>(std::floor(local_position[2]/ layer_widths[1]));

		return detID(module_index, layer_number, x_index, z_index);

	} //GetDetID

	//pointer Hit
	template<class Hit>
	detID GetDetID(Hit _hit){
		return GetDetID(_hit->x, _hit->y, _hit->z);
	}

	//for vector::Vector of position
	detID GetDetID(vector::Vector _hit){
		return GetDetID(_hit.x, _hit.y, _hit.z);
	}

	detID GetDetID(std::vector<double> _hit){
		return GetDetID(_hit[0], _hit[1], _hit[2]);
	}

	std::vector<double> GetCenterFloor(detID _id){
		return _floor.GetCenter(_id);
	}

    std::vector<double> GetCenterWall(detID _id){
        return _wall.GetCenter(_id);
    }

	std::vector<double> GetCenter(detID _id){
		if (_id.IsNull()) {
			std::cout << "detID is null" << std::endl;
			return {};
		}

		if (_id.isFloorElement) return GetCenterFloor(_id);
        if (_id.isWallElement) return GetCenterWall(_id);
        
		auto module = module_list[_id.moduleIndex];
		auto layer = layer_list[_id.layerIndex];

		std::vector<double> module_layer_center = {module->cx, layer->center, module->cz};
		std::vector<double> widths = layer->widths();


		module_layer_center[0] += widths[0]*static_cast<double>(_id.xIndex) + widths[0]/2.0;
		module_layer_center[2] += widths[1]*static_cast<double>(_id.zIndex) + widths[1]/2.0;


		if (GetDetID(module_layer_center[0], module_layer_center[1], module_layer_center[2]) == _id) return module_layer_center;
		else{
			std::cout << "WARNING: GEOMETRY MISMATCH! " << std::endl;
			std::cout << "Make sure GetDetID and GetCenter(detID) are consistent!  " << std::endl;

			_id.Print();
			std::cout << "module center: " << module->cx << ", " << module->cz << std::endl;
			std::cout << "widths : " << widths[0] << ", " << widths[1] << std::endl;
			std::cout << "calculated: " << module_layer_center[0] << ", " << module_layer_center[2] << std::endl;


			return module_layer_center;
		}
	}


	bool inBox(double x, double y, double z) {

		if (detector::x_min < x < detector::x_max) {
			if (detector::y_min < y < detector::y_max) {
				if (detector::z_min < z < detector::z_max) {
					return true;
				};
			};
		};

		return false;
	}



}; //Geometry Class




#endif
