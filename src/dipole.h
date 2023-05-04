/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2020 U. Lehnert
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
 * \class HardEdgeDipole
 * \brief A magnetic dipole
 * @author Ulf Lehnert
 * @date 10.7.2020
 * 
 * This lattice element describes a region of space with a constant dipole magnetic field.
 * The region is limited by two infinite planes each defined by
 * one point and the outside normal direction.
 * Field magnitude and direction can be defined with arbitrary values.
 * The reported field is constant within the region and drops to zero immediately outside.
 */
class HardEdgeDipole : public LocalizedField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet produce any field output.
     */
    HardEdgeDipole();

    /*! This constructor places the field at a given position in lab space.
     *  @param B field vector
     *  @param p1 position of the first limiting plane
     *  @param n1 outside normal of the first limiting plane
     *  @param p2 position of the second limiting plane
     *  @param n1 outside normal of the second limiting plane
     */
    HardEdgeDipole(Vector B, Vector p1, Vector n1, Vector p2, Vector n2);

    /*! This constructor takes the information from an XML node
     *  describing all field properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    HardEdgeDipole(const pugi::xml_node node, InputParser *parser);

    //! Get the field vector [T]
    Vector GetB() { return B_val; };

private:

    ElMagField LocalField(double t, Vector X) override;

    Vector  B_val;                              //! field value [T]
    Vector  p_1;                                //! position of the first limiting plane
    Vector  n_1;                                //! outside normal of the first limiting plane
    Vector  p_2;                                //! position of the second limiting plane
    Vector  n_2;                                //! outside normal of the second limiting plane

};


/*!
 * \class SoftEdgeDipole
 * \brief A magnetic dipole
 * @author Ulf Lehnert
 * @date 28.9.2021
 * 
 * This lattice element describes a region of space with a constant dipole magnetic field.
 * The region is limited by two infinite planes each defined by
 * one point and the outside normal direction.
 * Field magnitude and direction can be defined with arbitrary values.
 * The reported field is constant within the region and drops to zero outside.
 * The edges are modeled with an ArcTan shape, the given transition length
 * corresponds to the transition from 75% to 25% of the peak field.
 */
class SoftEdgeDipole : public LocalizedField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet produce any field output.
     */
    SoftEdgeDipole();

    /*! This constructor places the field at a given position in lab space.
     *  @param B field vector
     *  @param p1 position of the first limiting plane
     *  @param n1 outside normal of the first limiting plane
     *  @param t1 transition length at the first field boudary
     *  @param p2 position of the second limiting plane
     *  @param n2 outside normal of the second limiting plane
     *  @param t2 transition length at the second field boudary
     */
    SoftEdgeDipole(Vector B, Vector p1, Vector n1, double t1, Vector p2, Vector n2, double t2);

    /*! This constructor takes the information from an XML node
     *  describing all field properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    SoftEdgeDipole(const pugi::xml_node node, InputParser *parser);

    //! Get the field vector [T]
    Vector GetB() { return B_val; };

private:

    ElMagField LocalField(double t, Vector X) override;

    Vector  B_val;                              //! field value [T]
    Vector  p_1;                                //! position of the first limiting plane
    Vector  n_1;                                //! outside normal of the first limiting plane
    double  t_1;                                //! transition length at the first field boudary
    Vector  p_2;                                //! position of the second limiting plane
    Vector  n_2;                                //! outside normal of the second limiting plane
    double  t_2;                                //! transition length at the second field boudary

};


