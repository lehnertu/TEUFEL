#ifndef GENGRID_H
#define GENGRID_H
#include "vector.h"
#include <math.h>


using namespace std;

//class to generate a very simple planar grid centered around 0,0,z.

class GenGrid
{
public:
    //constructor
    GenGrid();

     //generates a grid parallel to x-y plane.
     //currently intended for seeing the flux at a screen placed parallel to x-y plane.
     //creates the grid extending from -x/2,-y/2 to +x/2,+y/2 and located at z
    void GenPlanarGrid(double x, double y, double z,int NumGridpoints);

    //Gets the position vector of the ith tile in the grid.
    // The origin point is at the beginning of undulator.
    Vector GetGridPoints(int i);

    //Gets the number of tiles in the grid generated.
    int GetGridSize();

    //Gets the area of a tile in the grid.
    double GetMeshArea();

   // Gets a normal vector to the grid -
   //which is presently the unit vector in z direction.
    Vector GetNormalVector();
private:
    Vector *r;
    int Numpoints;
    int TruePoints;
    double MeshArea;
    Vector MeshNormal;
};

#endif // GENGRID_H

