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

#include "vector.h"
#include "fields.h"

/*!
 * \class PlanarUndulator
 * \brief A planar undulator
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This planar undulator has a sinusoidal field with period \f$\lambda\f$ along the z-axis.
 * The main field component points in y direction with the peak value given
 * in the constructor. 
 * The field infinitely extends along the x-axis with constant value.
 * In y direction the field is also periodic with the same period as in z-axis but
 * the functional dependece is cosh() so there are poles at \f$\pm\lambda/2\f$.
 * 
 * In its own local coordinate system the undulator ist centered about the origin
 * extending \f$N\,\lambda/2\f$ in both directions along the z axis.
 * The edges are modeled with a linear ramp extending \f$\pm\lambda/2\f$ about
 * the z-axis values of the entrance and exit. This ensures an approximately
 * symmetric oszillation about the initial trajectory in x direction.
 * 
 * \todo The shift of origin is not yet handled.
 */
class PlanarUndulator : public ExternalField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet produce any field output.
     */
    PlanarUndulator();

    /*! This constructor places the undulator at a given position in lab space
     *  and initalizes all variables with sane values. This will not yet produce any field output.
     */
    PlanarUndulator(Vector pos);

    /* Specify the magnetic field.<br>
     * This is called with values B=0, lambda=1.0, N=1 by the constructors
     */
    void Setup(
	double B,                         // peak field [T]
        double lambda,                    // undulator period [m]
        int    N                          // number of undulator periods
        );

    double  GetBPeak();
    double  GetLambdaU();
    int     GetNPeriods();
    double  GetKpeak();
    double  GetKrms();

private:

    ElMagField LocalField(double t, Vector X);

    double  BPeak;                              // peak field [T]
    double  LambdaU;                            // undulator period [m]
    int     NPeriods;                           // number of undulator periods
    double  Krms;                               // undulator parameter
    double  ky, kz;                             // undulator periodicity [1/m]
};
