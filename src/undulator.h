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

#include "pugixml.hpp"
#include "vector.h"
#include "fields.h"
#include "parser.h"

/*!
 * \class PlanarUndulator
 * \brief A planar undulator
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This planar undulator has a sinusoidal field with period \f$\lambda\f$ along the z-axis.
 * The main field component points in y direction with the peak value given
 * in the setup. 
 * The field infinitely extends along the x-axis with constant value.
 * In y direction the functional dependece of the field is cosh() so the field diverges far from the axis.
 * 
 * In its own local coordinate system the undulator ist centered about the origin
 * extending \f$N\,\lambda/2\f$ in both directions along the z axis.
 * The edges are modeled with a linear ramp extending \f$\pm\lambda/2\f$ about
 * the z-axis values of the entrance and exit. This ensures an approximately
 * symmetric oszillation about the initial trajectory in x direction.
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

    /*! This constructor takes the information from an XML node
     *  describing all undulator properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    PlanarUndulator(const pugi::xml_node node, InputParser *parser);

    /* Specify the magnetic field.<br>
     * This is called with values B=0, lambda=1.0, N=1 by the constructors.
     */
    void Setup(
	double B,                          // peak field [T]
        double lambda,                 // undulator period [m]
        int    N                       // number of undulator periods
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


/*!
 * \class TransverseGradientUndulator
 * \brief A planar undulator with an additional transverse field gradient
 * @author Ulf Lehnert
 * @date 27.5.2019
 * 
 * This planar undulator has a sinusoidal field with period \f$\lambda\f$ along the z-axis.
 * The main field component points in y direction with the peak value on axis given in the setup. 
 * The field is supposed to change along the x direction with a given linear gradient.
 * To achieve this, a quadrupole-like field is added which
 * has a \f$\sin(k_x x)\f$ dependence on the horizontal coordinate and has zero divergence and curl.
 * The magnitude of this field is given by the scaling factor \f$k_x\f$.
 * When extrapolated linearly from the axis the field doubles at \f$x=1/k_x\f$ and reduces to zero at \f$x=-1/k_x\f$.
 * (\f$k_x\f$ can be negative as well.) The field oscillates when going far from the axis but
 * is approximately linear from 50% to 150% of the on-axis field.
 * 
 * In y direction the the functional dependece of the field is \f$\cosh(k_y y)\f$ so the field diverges far from the axis.
 * To satisfy the Maxwell equations \f$k_y\f$ is chosen as \f$k_y^2 = k_z^2 + k_x^2\f$.
 * 
 * @TODO To compensate for the beam deflection due to the gradient a dipole field
 * can be added. Along with the dipole comes a quadrupole which strength is
 * automatically adjusted such that the transverse gradient of the dipole strength
 * matches the transverse gradient of the undulator field. This is the expected field profile
 * when the iron poles of an undulator are additionally excited with a dipole coil.
 * 
 * In its own local coordinate system the undulator ist centered about the origin
 * extending \f$N\,\lambda/2\f$ in both directions along the z axis.
 * The edges are modeled with a linear ramp extending \f$\pm\lambda/2\f$ about
 * the z-axis values of the entrance and exit. This ensures an approximately
 * symmetric oszillation about the initial trajectory in x direction.
 */
class TransverseGradientUndulator : public ExternalField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet produce any field output.
     */
    TransverseGradientUndulator();

    /*! This constructor places the undulator at a given position in lab space
     *  and initalizes all variables with sane values. This will not yet produce any field output.
     */
    TransverseGradientUndulator(Vector pos);

    /*! This constructor takes the information from an XML node
     *  describing all undulator properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    TransverseGradientUndulator(const pugi::xml_node node, InputParser *parser);

    /* Specify the magnetic field.<br>
     * This is called with values B=0, lambda=1.0, N=1 by the constructors.
     */
    void Setup(
        double B,                       // peak field [T]
        double grad,                    // transverse gradient scaling factor kx [1/m]
        double modu,                    // pole-strength modulation
        double lambda,                  // undulator period [m]
        int    N                        // number of undulator periods
        );

    double  GetBPeak();
    double  GetLambdaU();
    int     GetNPeriods();
    double  GetKpeak();
    double  GetKrms();

private:

    ElMagField LocalField(double t, Vector X);

    double  BPeak;                              // peak field [T]
    double  Modulation;                         // pole strength (fractional) deviation of positive and negative poles
    double  LambdaU;                            // undulator period [m]
    int     NPeriods;                           // number of undulator periods
    double  Krms;                               // undulator parameter
    double  kx;                                 // transverse scaling factor [1/m]
    double  ky, kz;                             // undulator periodicity [1/m]
};
