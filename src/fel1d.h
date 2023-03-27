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
 * \class FEL_1D
 * \brief An electromagnetic field interacting with the beam
 * @author Ulf Lehnert
 * @date 27.3.2023
 * 
 * This is an electromagnetic field described on a 1D grid.
 * The grid propagates with the wave at the speed of light.
 * The transverse extent of the field is defined by a prescribed optical mode definition.
 * One polarization direction perpendicular to the direction of propagation is described.
 * 
 * The index 0 grid cell is placed at the given position in space and time.
 * The grid extends spatially opposite to the direction of propagation, so,
 * all grid cells subsequentially pass through the origin point.
 * The cell size must be the time-step divided by the speed of light for
 * the diffraction algorithm to work properly.
 */
class FEL_1D : public InteractionField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet produce any fields.
     */
    FEL_1D(
        double time_step,
        int number_steps
        );

    /*! This constructor takes the information from an XML node
     *  describing all field properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    FEL_1D(
        double time_step,
        const pugi::xml_node node,
        InputParser *parser
        );

    /*! All derived classes from GeneralField must provide a destructor */
    virtual ~FEL_1D();

    /*!
     * The electromagnetic field at a given time and point in space.
     * 
     * The coordinates [m] and the time [s] refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     */
    virtual ElMagField Field(double t, Vector X);

private:

    double  N_field;                            // number of the field grid cells
    double  dt;                                 // time step [s]

};

