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
 *    This class provides means to create and handle certain distributions
 *    of particles in phase space. It is meant for the creation of 6D
 *    coordinate distributions which can be used by the bunch class
 *    to initialize particles for tracking.
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
    
    /*! Join in another distribution.
     *  This only works if both have the same number of dimensions.
     */
    void join(const Distribution *other);
    
    //! report number of dimensions
    int getDIM();
    
    //! report number of points
    int getNOP();
    
    /*! generate a gaussian distribution for one of the dimensions
     *  @todo: This should generate a strictly sequential Gaussian.
     *  The present code shoud be renamed generateNormalDist()
     */
    void generateGaussian(double mean, double sigma, int dim);
    
    /*! Add a correlation between two axis.
     *  The value of the independent coordinate multiplied with a factor
     *  is added to the dependent axis coordinate value.
     */
    void addCorrelation(int independent, int dependent, double factor);
    
    /*! Get one coordinate of one particle with given index.
     *  Out of range indices will not lead to errors, just return zero
     */
    double getCoordinate(int index, int dim);
    
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
    
    Creating a bunch does not define inital coordinates of the particles.
    This is done when initalizing the tracking algorithm (e.g. InitVay() ).
 */
class Bunch
{

public:

    /*! default constructor: creates an empty bunch .
    */	
    Bunch();
    
    /*! create a bunch of given number of particles each having
     * given charge and mass.
     * 
     * The particles have no trajectory data, therefore, also no initial coordinates
     */
    Bunch(int N, double charge, double mass);
    
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

    /*! @brief Setup for tracking the whole bunch using the Vay algorithm.
     * 
     * The inital coordinates of all particles are defined by the given distribution.
     * 
     * See ChargedParticle::InitVay for details
     * 
     * The given distribution contains 3 position coordinates [m],
     * 3 momentum coordinates [] and the start tim [s].
     * If the given distribution has less than 7 dimensions the missing
     * coordinates are initalized as 0 which can be used to start all
     * particles at zero time.
     * 
     * @todo Tracking particles which do not start all at the same time
     * is not yet supported.
     * 
     * @param dist the inital coordinates of all particles
     * @param tstep the length of the time step.
     * @param field the field through which the particle will be tracked
     */
    void InitVay(Distribution *dist,
		 double tstep,
		 GeneralField *field);

    /*! @brief Perform one tracking step using the Vay algorithm.
     * 
     * See ChargedParticle::InitVay for details
     * 
     * @param field the field through which the particle is tracked
     */
    void StepVay(GeneralField* field);

    /*! Dump all particle information belonging to a given observation time
     *  into an SDDS file.
     * 
     * returns values for error checks:
     *	 
     *	0  -  successfully Written the file\n
     *	1  -  error in SDDS_InitializeOutput \n
     *	2  -  error in SDDS_DefineSimpleParameter \n
     *	3  -  error in SDDS_DefineColumn \n
     *	4  -  error in SDDS_WriteLayout \n
     *	5  -  error in SDDS_StartPage \n
     *	6  -  error in SDDS_SetParameters \n
     *	7  -  error in SDDS_SetRowValues \n
     *	8  -  error in SDDS_WritePage \n
     *	9  -  error in SDDS_Terminate \n
     * 
     */
    int WriteWatchPointSDDS(double time,
			    const char *filename);

    /*! Dump all particle information belonging to a given observation time
     *  into an HDF5 file.
     */
    int WriteWatchPointHDF5(double time,
			    const char *filename);


private:

    //! Number of Particles in the bunch
    int NOP;

    //! time step for tracking - this will remain constant after being set at the start of tracking
    double dt;

    //! we store references to all particles
    vector<ChargedParticle*> P;

};