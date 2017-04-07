#ifndef GENGRID_H
#define GENGRID_H
#include "vector.h"
#include <math.h>
#include "particle.h"
#include "SDDS.h"
#include "bunch.h"
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

    /*!

	 write SDDS file to output file named "radiation@grid.sdds

	 writes data for every time step in different page

	 Arguments have the following meaning:\n
	 Vector Robs -> Observation Point of the Radiation\n
	 time_begin  -> starting point or expected arrival time of the signal\n
	 time_end    -> expected time at which the last wavefront would leave\n
	 NumberOfPoints -> Total Number of Sample Points for the radiation\n	

	 use sddsquery trajectory.sdds to view the format of dataset
	 
	 routine returns values with following meaning:
	 
	0  ->  Successfully Written the file \n
	1  ->  Error in Initializing the Output dataset \n
	2  ->  Error in Defining the parameters \n
	3  ->  Error in Defining Columns describing the data to be followed \n
	4  ->  Error in Writing the layout of the data structure in the sdds file \n
	5  ->  Error in Starting a New Page of the file  \n
	5  ->  Error in setting the values of the parameters \n
	6  ->  Error in setting the row values i.e. the data belonging to column \n
	7  ->  Error in Writing the page that was successfully initialized \n
	8  ->  Error in Terminating the data flow to the sdds file \n

	*/
    int SDDSRadiationAtGrid(Bunch *bunch, double time_begin, double time_end, int NumberOfPoints );



private:
    Vector *r;
    int Numpoints;
    int TruePoints;
    double MeshArea;
    Vector MeshNormal;
};

#endif // GENGRID_H

