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

#include <complex>

#include "pugixml.hpp"
#include "vector.h"
#include "fields.h"
#include "parser.h"

/*!
 * \class GaussianWave
 * \brief A Gaussian electroagnetic wave
 * @author Ulf Lehnert
 * @date 9.4.2019
 * 
 * This class describes a gaussian mode of electromagnetic radiation
 * describing a propagating paraxial wave. The origin of the field is
 * the mode waist. The propagation direction is along the z-axis and
 * the polarization of the electic field is along the x-axis.
 * Possitions and directions can be changed by transformations
 * given in the constructors applied through the ExternalField base class.
 * The transverse mode seize is defined through the Rayleigh length.
 * 
 * The radiation power is constant all along the beam with the
 * field amplitude varying accordin to the mode size. It is defined through
 * the amlitude at the waist given as a complex value of the electric field  strength.
 * The real part describes the \f$ cos(\omega t) \f$ terms, the imaginary part
 * the  \f$ sin(\omega t) \f$ terms, thus, defining the phase at the waist.
 * 
 * The total power is \f$ P = \frac{\pi}{2} I_0 \omega_0^2 \f$
 * with \f$ I_0 = |E_0|^2 / (2 \eta) \f$ the peak intensity on axis.
 * \f$ \eta = 1 / (\epsilon_0 c) \f$ is the vacuum impedance.
 */
class GaussianWave : public ExternalField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet
     *  produce any field output as the amplitude is zero.
     */
    GaussianWave();

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet
     *  produce any field output as the amplitude is zero.
     */
    GaussianWave(Vector pos);

    /*! This constructor takes the information from an XML node
     *  describing all wave properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    GaussianWave(const pugi::xml_node node, InputParser *parser);

    /*! Specify the electromagnetic field properties.<br>
     *  This is called by the constructors with following default values :
     *  - omega = 1.0
     *  - A = (0.0, 0.0)
     *  - range = 1.0
     */
    void Setup(
        double lambda,					//! wavelength [m]
        complex<double> amplitude,		//! the field amplitude at the waist center [V/m]
        double range					//! the Rayleigh length [m]
        );

	/*! Report the wavelength lambda [m]
	 */
	double getWavelength() { return lambda; }
	
private:

	/*! compute the electromagnetic fields in local coordinates */
    ElMagField LocalField(double t, Vector X);

    double lambda;						//! wavelength [m]
    double k;							//! wave number [1/m]
    double omega;						//! circular frequency [Hz]
    complex<double> A;					//! the field amplitude at the waist center [V/m]
    double zR;							//! the Rayleigh length [m]
    double w0;							//! mode size at waist [m]

};

