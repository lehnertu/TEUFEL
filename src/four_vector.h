/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2017 U. Lehnert
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#ifndef FOUR_VECTOR_H
#define FOUR_VECTOR_H
#include "vector.h"
#include "vector.cpp"
using namespace std;

class Four_Vector
{
	public:
		Four_Vector();
		Four_Vector(double t0, double x0, double y0, double z0);
		//for the four vectors, the first element has to be ct
		// t1 will be actual time
		// SpeedOfLight is multiplied to t0
		// Note the following:
		// Four Position = (ct,x,y,z)
		// four Momentum = (E,pxc,pyc,pzc);
		// special type of four vector
		// arguments: double E (Energy) Vector(p)(not multiplied by c)
		Four_Vector Four_Momentum(double E, Vector p);
		double t;
		double x;
		double y;
		double z;
		//Lorentz Transform from the lab frame to some frame
		//the frame has gamma beta equal to GB;
		//restricting the axes of frames to be parallel
		Four_Vector LorentzTransform(Four_Vector X, Vector GB);
		Four_Vector InverseLorentzTransform();
		
		// sum and difference of four two vectors
    		Four_Vector operator+ (Four_Vector);
    		Four_Vector& operator+= (Four_Vector);
    		Four_Vector operator- (Four_Vector);

    		// negative vector
    		Four_Vector operator- ();

   		 // multiplication of a vector with a factor
   		Four_Vector operator* (double faktor);
    		Four_Vector& operator*= (double faktor);
    		Four_Vector operator/ (double faktor);
    		Four_Vector& operator/= (double faktor);

    		// absolute value of a vector
    		double norm();

    		// squared absolute value of a vector
    		double abs2nd();

    		// make a vector unit length
    		void normalize();

		//transform a four vector X into frame S1
		Four_Vector LorentzTransform(Vector S1);

		//inverse transform vector X into frame S1
		Four_Vector InverseLorentzTransform(Vector S1);


		//find the dot product of twoo vectors
		double dot(Four_Vector a, Four_Vector b);

	private:
		//LorentZ Transformation Matrix
		// takes the gamma beta of the moving frame
		// creates a 4 by 4 matrix of the transformation elements
		// Vector S1 = Vector(gamma*betax, gamma*betay,gamma*betaz)
		double TransferMatrix(int i, int j,Vector S1);
		double InverseTransferMatrix(int i, int j,Vector S1);
		
		


};

#endif
