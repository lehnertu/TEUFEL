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

#include "fields.h"
#include "bunch.h"
#include "vector.h"
#include <vector>

using namespace std;

/*!
    \class Beam
    \brief Ensemble of bunches
 
    @author Ulf Lehnert
    @date 19.9.2017
    
    This is a container holding a number of bunches.
    The main purpose is to parallelize the tracking of several buches but
    also tracking of several distinct fractions of a bunch is possible.
    The interface almost exactly resembles that of a single bunch.
 */
class Beam
{

public:

    /*! default constructor: creates an empty Beam.
     */	
    Beam();
    
    /*!
     * Destructor:
     * Deleting the bunch also deletes all contained bunches.
     * This in turn deletes all particles belonging to the bunches.
     * So, normally the beam is the only object that needs to be deleted.
     */
    ~Beam();
    
    /*!
     * Add a bunch to the beam.
     * The beam then "owns" the bunch and will delete it when destructed.
     */
    void Add(Bunch *bunch);
    
    //! Report the number of bunches in the beam.
    int getNOB();	
    
    //! Report the total charge of the particles contained in the beam.
    double getTotalCharge();
    
    /*! @brief Setup for tracking the whole beam using the Vay algorithm.
     * 
     * See ChargedParticle::InitVay for details
     * 
     * @todo Tracking particles which do not start all at the same time
     * is not yet supported.
     * 
     * @param tstep the length of the time step.
     * @param field the field through which the beam will be tracked
     */
    void InitVay(double tstep,
		 GeneralField *field);
    
    /*! @brief Perform one tracking step using the Vay algorithm.
     * 
     * See ChargedParticle::StepVay for details
     * 
     * @param field the field through which the beam is tracked
     */
    void StepVay(GeneralField *field);
    
    /*! Dump all particle information into an HDF5 file.
     *  The written quantities include:
     *  - time t
     *  - position x,y,z
     *  - momentum px,py,pz,p (beta*gamma)
     *  - angles xp,yp (px/pz, py/pz)
     *  - gamma
     */
    int WriteWatchPointHDF5(const char *filename);
    
private:

    //! Number of bunches in the beam.
    int NOB;

    //! we store references (pointers) to all bunches
    vector<Bunch*> B;

    /*! Time step for tracking - this will remain constant
     * after being set at the start of tracking.
     * @todo do we actually need it?
     */
    double dt;

};
