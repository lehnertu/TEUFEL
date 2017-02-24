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

#ifndef VECTOR_H
#define VECTOR_H

using namespace std;

class Vector
{

  public:

    Vector();
    Vector(double x0, double y0, double z0);

    double x;
    double y;
    double z;

    // sum and difference of two vectors
    Vector operator+ (Vector v2);
    Vector& operator+= (Vector v2);
    Vector operator- (Vector v2);

    // negative vector
    Vector operator- ();

    // multiplication of a vector with a factor
    Vector operator* (double faktor);
    Vector& operator*= (double faktor);
    Vector operator/ (double faktor);
    Vector& operator/= (double faktor);

    // absolute value of a vector
    double norm();

    // squared absolute value of a vector
    double abs2nd();

    // make a vector unit length
    void normalize();

  private:

};

// dot product of two vectors
double dot(Vector a, Vector b);

// cross product of two vectors
Vector cross(Vector a, Vector b);

#endif
