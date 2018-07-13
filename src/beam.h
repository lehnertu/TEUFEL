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
    
    //! Report the number of bunches in the beam.
    int getNOB();       

    /*!
     * Get a pointer to a particular bunch.
     * Returns a NULL pointer if the index is out of range.
     */
    Bunch* getBunch(int i);
    
    //! Report the total number of particles contained in the beam.
    int getNOP();
    
    //! Report the total charge of the particles contained in the beam.
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
    
    /*! Compute the electromagnetic field radiated by the bunch
     * seen at the observation point. The field is given in time domain
     * starting at t0 with NOTS equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * So the center (reference) time for the first sample is \f$t0+dt/2\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * 
     * This method sums up the field generated by the individual particles.
     * All field contributions are just added to the field referenced
     * by ObservationField. The caller may choose to start with a zero field
     * or to use this to sum up the fields of multiple time steps
     * or multiple bunches.
     * 
     * This method linearly interpolates over time, assuming the time
     * step for the trajectory integration is choosen small enough.
     * Thus, the field-time integral is exactly preserved.
     * 
     * \param[in] ObservationPoint The position [m] of the observer.
     * \param[in] t0 Start time [s] of the trace (at the observer).
     * \param[in] dt Duration of one time segment [s] of the observation trace.
     * \param[in] nots Number of time segments.
     * \param[out] ObservationField Field integrals at the observation point.
     */
    void integrateFieldTrace(
	Vector ObservationPoint,
	double t0,
	double dt,
	int nots,
	std::vector<ElMagField> *ObservationField);
    
private:

    //! Number of bunches in the beam.
    int NOB;

    //! we store references (pointers) to all bunches
    vector<Bunch*> B;

    /*! Time step for tracking - this will remain constant
     * after being set at the start of tracking.
     */
    double dt;

    //! numer of time steps
    int NOTS;
    
    //! Method used for tracking
    TrackingMethod tracker;
    
};
