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
#include "particle.h"
#include "vector.h"
#include <vector>

using namespace std;

/*!
 *    \class Distribution
 *    \brief Generate arrays of particle coordinates.
 * 
 *    @author Ulf Lehnert
 *    @date 17.5.2017
 *    
 *    This class provides means to create certain distributions of particles
 *    in phase space.
 */
class Distribution
{

public:
    
    /*! By default the distribution is initalized with particles
     *  uniformly distributed over the range [0.0, 1.0) in all coordinates.
     *  The number of dimensions of the phase space and the number of particles
     *  to be created must be given.
     */
    Distribution(int dim, int nop);
    
    //! Destructor only used to free the memory
    ~Distribution();
    
private:
    
    //! the number of dimensions
    int DIM;
    
    //! the number of particles
    int NOP;
    
    //! the array of coordinates
    double* A;
    
};

    /*!
    \class Bunch
    \brief Ensemble of particles
 
    @author Ulf Lehnert, Vipul Joshi
    @date 10.5.2017
    
    This is a container holding a number of particles. Tracking particles
    and computation of radiated fields are provided for all particles together.
 */
class Bunch
{

public:

    /*! default constructor: creates an empty bunch .
    */	
    Bunch();

    /*!
     * copy constructor:
     * Create a copy of an existing bunch, thereby, creating copies of all particles
     * of the original bunch. The original bunch is not altered.
     */
    Bunch(Bunch* origin);

    /*!
     * Destructor:
     * Deleting the bunch also deletes all contained particles
     */
    ~Bunch();

    /*!
     * Add a particle to the bunch.
     * The bunch then "owns" the particle and will delete it when destructed.
     */
    void Add(ChargedParticle *part);

    //! Report the number of particles in the bunch.
    int getNOP();	

    //! Report the total charge of the particles contained in the bunch.
    double getTotalCharge();

    //! Get a pointer to a particle from its index within the bunch
    ChargedParticle* getParticle(int i);

private:

    //! Number of Particles in the bunch
    int NOP;

    //! we store references to all particles
    vector<ChargedParticle*> P;

};