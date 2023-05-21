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

#include "beam.fwd.h"
#include "fields.fwd.h"
#include "global.h"
#include "bunch.h"
#include "vector.h"
#include <vector>

using namespace std;

//! choice of the tracking method
typedef enum keyword { TRACKING_NONE, TRACKING_EULER, TRACKING_VAY } TrackingMethod;

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
     * Deleting the beam also deletes all contained bunches.
     * This in turn deletes all particles belonging to the bunches.
     * So, normally the beam is the only object that needs to be deleted.
     */
    ~Beam();
    
    /*!
     * Add a bunch to the beam.
     * The beam then "owns" the bunch and will delete it when destructed.
     */
    void Add(Bunch *bunch);
    
    /*! Remove all bunches and particles belonging to this beam.
     */
    void clear();
    
    /*! Remove all trajectory information from particles
     *  belonging to this beam
     */
    void clearTrajectories();
    
    /*! For all particles belonging to this beam pre-allocate a number of trajectory points.
     *  Note - this must be one more than the number of tracking steps
     */
    void preAllocate(int nTraj);
    
    //! Report the number of bunches in the beam.
    int getNOB();       

    /*!
     * Get a pointer to a particular bunch.
     * Returns a NULL pointer if the index is out of range.
     */
    Bunch* getBunch(int i);
    
    //! Report the total number of particles contained in the beam.
    int getNOP();
    
    /*! Report the total charge of the particles contained in the beam
     *  in units of elementary charge (i.e. simulated number of electrons)
     */
    double getTotalCharge();
    
    //! Report the used tracking method
    TrackingMethod getTrackingMethod() { return tracker; };
    
    //! Set the tracking method to be used.
    void setTrackingMethod(TrackingMethod method) { tracker = method; };
    
    //! Report the used tracking time step
    double getTimeStep() { return dt; };
    
    //! Set the tracking time step to be used.
    void setTimeStep(double step) { dt = step; };
    
    //! Report the number of time steps
    int getNOTS() { return NOTS; };
    
    //! Set the number of time steps to be used for tracking.
    void setNOTS(int n) { NOTS = n; };

    /*!
     *  Generic setup method to initlize the tracking of all particles in the beam.
     *  Depending on the choosen stepper method the appropriate init methods
     *  are called.
     * 
     *  @param field the field through which the beam will be tracked
     * 
     *  In case of failing sanity checks an IOexception is thrown.
     */
    void setupTracking(GeneralField *field);
    
    /*!
     *  Generic setup method to to do oe tracking step of all particles in the beam.
     *  Depending on the choosen stepper method the appropriate stepper is called.
     * 
     *  @param field the field through which the beam is tracked
     * 
     *  In case of failing sanity checks an IOexception is thrown.
     */
    void doStep(GeneralField *field);
    
    /*! @brief Compute the buffer size for one step.
     * 
     * @return the number of doubles needed
     */
    int getStepBufferSize();
    
    /*! @brief Buffer particle coordinates.
     * 
     * Particle coordinates (time,position,momentum,acceleration)
     * of all particles belonging to the beam are stored into one buffer.
     * The necessary buffer size is given by getStepBufferSize().
     * There is no memory check for the buffer size performed.
     * 
     * @param[in] buffer memory pointer of the buffer
     */
    void bufferStep(double *buffer);
    
    /*! @brief Set particle coordinates from buffer.
     * 
     * Particle coordinates (time,position,momentum,acceleration)
     * of all particles belonging to the beam are read from one buffer.
     * The necessary buffer size is given by getStepBufferSize().
     * There is no memory check for the buffer size performed.
     * 
     * @param[in] buffer memory pointer of the buffer
     */
    void setStepFromBuffer(double *buffer);
    
    /*! @brief Compute the buffer size for one set of probe particles.
     *
     *  @param[in] the number of particles requested
     *  @return the number of doubles needed
     */
    unsigned int getProbeBufferSize(unsigned int number);

    /*! @brief Buffer particle information.
     * 
     *  Particle coordinates (time,position,momentum,acceleration)
     *  perceived fields (E,B)
     *
     *  The necessary buffer size is given by getProbeBufferSize().
     *  There is no memory check for the buffer size performed.
     * 
     *  @param[in] buffer memory pointer of the buffer
     *  @param[in] the number of particles requested
     */
    void bufferProbe(double *buffer, unsigned int number);

    /*! Dump all particle information into an HDF5 file.
     *  The written quantities include:
     *  - time t
     *  - position x,y,z
     *  - momentum px,py,pz,p (beta*gamma)
     *  - angles xp,yp (px/pz, py/pz)
     *  - gamma
     * 
     *  @return number of particles written, -1 in case of an error
     */
    int WriteWatchPointHDF5(std::string filename);
    
    /*! Dump all particle information into an SDDS file.
     *  The particles are transported to the crossing with an average plane perpendicular to the direction of beam propagation.
     *  The SDDS file uses a left-handed cartesian coordinate system x-y-s (s is the propagation direction).
     *
     *  A number of parameters is written to the header:
     *  - Particles - the total number of particles in the file
     *  - Charge - the total charge
     *  - pCentral - average beta*gamma (reference momentum)
     *  - sdir_x, sdir_y, sdir_z - the propagation direction of the beam in laboratory coordinates
     *  - xdir_x, xdir_y, xdir_z - the x direction of the file in laboratory coordinates
     *  - ydir_x, ydir_y, ydir_z - the y direction of the file in laboratory coordinates
     *
     *  The written columns include:
     *  - time t
     *  - position x,y
     *  - momentum p (exactly beta*gamma = p/m_e c)
     *  - angles xp,yp (dx/ds, dy/ds)
     * 
     * returns values for error checks:
     *	 
     *	0  -  file successfully written\n
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
    int WriteWatchPointSDDS(std::string filename);

    /*! @brief Get the field generated by a beam of particles.
     *
     *  Compute the electromagnetic field radiated by the particles of a beam
     *  as seen at the observation point at the given time.
     *  This method sums up the field generated by the individual particles.
     *  
     *  @param[in] obs_pos the position of the observer
     *  @param[in] obs_time the time of observation
     *
     *  @return the observed electromagnetic field
     *
     *  @throws RANGEexception
     *  An exception is thrown if any particle is not properly initialized.
     *  <br>
     */
    ElMagField RetardedField(double obs_time, Vector obs_pos);

    /*! Compute the electromagnetic field radiated by the bunch
     * seen at the observation point. The field is given in time domain
     * as defined in the FieldTrace object.
     * 
     * This method sums up the field generated by the individual particles.
     * All field contributions are just added to the trace
     * The caller may choose to start with a zero field
     * or to use this to sum up the fields of multiple bunches.
     * 
     * \param[in] ObservationPoint The position [m] of the observer.
     * \param[out] trace Electromagnetic field integrals at the observation point.
     */
    void integrateFieldTrace(
        Vector ObservationPoint,
        FieldTrace *trace);

private:

    /*! @brief Setup for tracking the whole beam using the Vay algorithm.
     * 
     * See ChargedParticle::InitVay for details
     * 
     * @todo Tracking particles which do not start all at the same time
     * is not yet supported.
     * 
     * @param field the field through which the beam will be tracked
     */
    void InitVay(GeneralField *field);
    
    /*! @brief Perform one tracking step using the Vay algorithm.
     * 
     * See ChargedParticle::StepVay for details
     * 
     * @param field the field through which the beam is tracked
     */
    void StepVay(GeneralField *field);
    
private:

    //! Number of bunches in the beam.
    int NOB;

    //! we store references (pointers) to all bunches
    vector<Bunch*> B;

    /*! Time step for tracking - this will remain constant
     * after being set at the start of tracking.
     */
    double dt;

    /*! number of time steps
     *  This is the number of time steps that have actually been performed during tracking the beam.
     *  The value is initially set to the requested number when parsing, but zeroed when
     *  tracking is initialized. Then it is incremented with every tracking step performed.
     *  Note - it is by 1 smaller than the number of trajectory points stored for every particle.
     */
    int NOTS;
    
    //! Method used for tracking
    TrackingMethod tracker;
    
};
