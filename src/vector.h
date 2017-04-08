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

#pragma once

using namespace std;

/*!
    \class Vector
    \brief 3D-Vector class
 
    @author Ulf Lehnert
    @date 10.2.2017

    \todo make components private
 */
class Vector
{

  public:

    /*!
	Default constructor: <br>
	All components are initalized as zero.
    */
    Vector();
    
    /*!
	Standard constructor: <br>
	All components are initalized with given values.
    */
    Vector(double x0, double y0, double z0);

    double x;
    double y;
    double z;

    /*! Sum of two vectors */
    Vector operator+ (Vector v2);

    /*! Accumulating sum of two vectors */
    Vector& operator+= (Vector v2);

    /*! Difference of two vectors */
    Vector operator- (Vector v2);

    /*! The negative of a vector */
    Vector operator- ();

    /*! Multiplication of the vector with a real factor */
    Vector operator* (double factor);

    /*! In-place multiplication of the vector with a real factor */
    Vector& operator*= (double factor);

    /*! Division of the vector by a real factor */
    Vector operator/ (double factor);

    /*! In-place division of the vector by a real factor */
    Vector& operator/= (double factor);

    /*! The absolute value (cartesian norm) of the vector */
    double norm();

    /*! The squared absolute value of the vector */
    double abs2nd();

    /*! Make the vector unit length */
    void normalize();

  private:

};

/*! The dot product of two vectors */
double dot(Vector a, Vector b);

/*! The cross product of two vectors */
Vector cross(Vector a, Vector b);
