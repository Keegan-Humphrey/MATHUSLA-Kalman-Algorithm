#include <cstdlib>
#include <iostream>
#include <iomanip>

#ifndef LIN_ALG_HH
#define LIN_ALG_HH


namespace vector{

class Vector {
public:
        double x, y, z;
        // 3D Coordinates of the Vector


        Vector(double x, double y, double z)
        {
                // Constructor
                this->x = x;
                this->y = y;
                this->z = z;
        }
        double Magnitude();

        //diagonal metric represented as vector, each entry scales the corresponding component [used for residuals]
        double Magnitude(Vector DiagonalMetric);

        Vector operator+(Vector v); // ADD 2 Vectors
        Vector operator-(Vector v); // Subtraction
        double operator^(Vector v); // Dot Product
        Vector operator*(Vector v); // Cross Product
        Vector Scale(double c); //Scalar Multiply
        std::vector<double> std(); //convert to std::vector<double>
    
        friend std::ostream& operator<<(std::ostream& out, const Vector& v);
        // To output the Vector
};




}; // namespace vector


#endif