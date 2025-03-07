/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2019 U. Lehnert
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "global.h"
#include "fields.h"
#include "vector.h"
#include "bunch.h"
#include "beam.h"

/*!
 * \class PlatePEC
 * \brief A perfectly conducting plate scattering the fields generated by tracked particles.
 * @author Ulf Lehnert
 * @date 14.12.2021
 * 
 * This class is intended to simulate metallic surfaces scattering the electromagnetic fields
 * radiated by the tracked particles. It has rectangular shape sub-divided in a grid.
 * The radiation impinging from tracked particles and other lattice elements is
 * recorded on-the-fly during tracking. Reflected fields are computed using
 * the Kirchhoff formula for arbitrary points in space and time.
 * These fields are considered lattice fields for tracking the particles.
 * This requires the field storage including the derivatives to be updated
 * every tracking step. Finally, the field traces can be stored in a screen-like file.
 */
template <class objectT>
class PlatePEC : public GeneralField
{

public:

    /*! Standard constructor:<br>
     * 
     * This class is templated and spezialized for
     * Bunch() or Beam() being observed.
     * 
     * \param[in] obj The object (bunch or beam) to be observed.
     * \param[in] filename The name of the file to write.
     */
    PlatePEC(objectT *obj, const char *filename);

    /*!
     * Destructor
     */
    virtual ~PlatePEC();
    
    /*!
     * The electromagnetic field at a given time and point in space.
     * 
     * The coordinates [m] and the time [s] refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     */
    ElMagField Field(double t, Vector X) override;

    /*! The source has advanced one time step.
     *  Compute and store the quantities of interest.
     */
    void update();

private:

    //! the observed beam object
    objectT *Source;
    
    //! the file name
    const char* fn;

    //! the number of stored datasets
    unsigned int NOTS;

};

