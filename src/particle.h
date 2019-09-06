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
#include "vector.h"

//! size of the buffer for particle serialization
#define PARTICLE_SERIALIZE_BUFSIZE 12

/*!
 * \class ChargedParticle
 * \brief Particle with electric charge.
 * @author Ulf Lehnert
 * @date 5.10.2017
 * 
 * This class holds the data of a single charged particle with given mass and charge.
 * For tracking this class provides a number of algorithms to perform a single time-step.
 * The trajectory is stored comprising time [s], position [m], momentum [] and acceleration [1/s] in
 * the laboratory frame.
 * 
 * All charged particles are surrounded by electromagnetic fields and
 * accelerated particle produce radiation. The fields produced by the particle
 * at a given position in laboratory frame can be computed
 * taking into account the full retardation effect.
 * 
 * @todo One could implement extrapolation of the trajectories with zero acceleration.
 * When implemented for negative times, this would allow to compute realistic space charge
 * fields from the beginning of tracking.
 */
class ChargedParticle
{
    
public:
    
    /*! Default constructor:<br>
     * Charge is set to electron charge, mass to electron mass.
     * No trajectory is allocated.
     */
    ChargedParticle();
    
    /*! Standard constructor:<br>
     * Charge is is given in units of elementary charge, mass im units of electron mass.
     * No trajectory is allocated.
     */
    ChargedParticle(double charge, double mass);
    
    /*! Copy constructor:<br>
     * Create a replication of a given particle.
     */
    ChargedParticle(const ChargedParticle *origin);
    
    /*! Constructor from buffer:<br>
     *  Take the data stored with serialize() and recreate this
     *  as a new particle including initTrajectory().
     *  There is no check of buffer size.
     */
    ChargedParticle(double *buffer);
    
    /*! Destructor: free all memory. */
    virtual ~ChargedParticle();
    
    /*! Return the number of points in the trajectory */
    int getNP() { return NP; };

    /*! Return the charge of the particle */
    double getCharge() { return Charge; };
    
    /*! Return the time of the last stored trajectory point
     *  or zero if there is no trajectory.
     */
    double getTime();
    
    /*! Return the position of the last stored trajectory point
     *  or zero if there is no trajectory.
     */
    Vector getPosition();
    
    /*! Return the momentum of the last stored trajectory point
     *  or zero if there is no trajectory.
     */
    Vector getMomentum();
    
    /*! Copy all current information about the particle into one buffer
     *  of double type. This can be used to transfer the particle to a different
     *  MPI node but looses the trajectory information.
     *  All stored information belongs to the last known trajectory point.
     *  The serialized properties include:
     *  @item double Charge
     *  @item double Mass
     *  @item double Time[0]
     *  @item double[3] X[0]
     *  @item double[3] P[0]
     *  @item double[3] A[0]
     *  The buffer must have a size of PARTICLE_SERIALIZE_BUFSIZE doubles.
     *  There is no check for sufficient buffer size.
     */
    void serialize(double *buffer);
    
    /*! Return the time [s] for one point of the trajectory */
    double TrajTime(int step);
    
    /*! Return the position [m] for one point of the trajectory */
    Vector TrajPoint(int step);
    
    /*! Return the dimensionless momentum \f$\beta*\gamma\f$ for one point of the trajectory */
    Vector TrajMomentum(int step);
    
    /*! Return the acceleration [1/s] \f$a = \frac{d}{dt}(\beta\gamma)\f$ for one point of the trajectory */
    Vector TrajAccel(int step);
    
    /*! Initialize a trajectory by defining its first point.
     *  Any possibly existing trajectory data is deleted.
     */
    void initTrajectory(double t, Vector x, Vector p, Vector a);
    
    /*! Mirror the trajectory about a given plane.
     *  The plane is defined by a vector origin and its normal vector.
     *  The charge is scaled by an additional given factor which is also applied to qm and qmt2.
     */
    void Mirror(Vector origin, Vector normal, double f_charge);
    
    /*! @brief Setup for tracking the particle using the Vay algorithm.
     * 
     *  The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     *  The initial position, momentum and time of the particle must have been
     *  setup before. An exception is thrown otherwise.
     * 
     *  The tracking is started from the last defined point of the existing
     *  trajectory. This allows to change the step size of tracking.
     *  The acceleration of the last trajectory point is overridden
     *  by the one computed from the particles motion inside the given field.
     * 
     *  @param tstep the length of the time step.
     *  @param field the field through which the particle is tracked
     *
     *  @throws RANGEexception
     *  <br>
     */
    void InitVay(double tstep,
		 GeneralField* field);

    /*! @brief Perform one tracking step using the Vay algorithm.
     * 
     *  The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     *  @param field the field through which the particle is tracked
     */
    void StepVay(GeneralField* field);
    
    /*! Particle phase space position at a given time.
     *  Interpolates between time steps after tracking.
     *  Can also be used when the tracking is only initalized to
     *  retreive the inital particle position at t=0.
     *
     *  @todo This should be extended with an extrapolation to negative times.
     */
    void CoordinatesAtTime(double time, Vector *position, Vector *momentum);
    
    /*! @brief Get the field generated by a particle.
     *
     *  Compute the electromagnetic field radiated by the particle
     *  as seen at the observation point at the given time.
     *  
     *  Given an observation point in space and time
     *  the one (source) point on the particle trajectory is computed
     *  from which a speed-of-light signal arrives just at the right time
     *  at the observation point. This is done by bisecting the
     *  trajectory until one step is found where the source point
     *  lies within. Within this step a linear interpolation is done.
     *
     *  If necessary the trajectory is extrapolated to negative times
     *  assuming the acceleration of the particle is zero at negative times.
     *  In that case the source point can be computed analytically.
     *
     *  Then the Lienard-Wiechert formula is used to compute the observed fields.
     *  At positive times if no signal from the known trajectory has reached
     *  the observation point zero field is returned.
     *
     *  @param[in] obs_pos the position of the observer
     *  @param[in] obs_time the time of observation
     *
     *  @return the observed electromagnetic field
     *
     *  @throws RANGEexception
     *  An exception is thrown if the particle is not properly initialized.
     *  <br>
     */
    ElMagField RetardedField(double obs_time, Vector obs_pos);
        
    /*! Compute the electromagnetic field radiated by the particle
     * seen at the observation point. The field is given in time domain
     * starting at t0 with NOTS equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * So the center (reference) time for the first sample is \f$t0+dt/2\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * 
     * This method sums up the field generated by the particle over all time steps.
     * Fields emitted before the start of the trajectory are not included.
     * The methods RetardedTime() and RetardedField() are called
     * to sequentially process all time steps stored in the trajectory.
     * All field contributions are just added to the field referenced
     * by ObservationField. The caller may choose to start with a zero field
     * or to use this to sum up the fields of multiple particles.
     * 
     * This method linearly interpolates over time, assuming the time
     * step for the trajectory integration is choosen small enough.
     * Thus, the field-time integral is exactly preserved when mapping
     * the retarded time-steps of the trajectory onto the observation trace.
     * 
     * \param[in] ObservationPoint The position [m] of the observer.
     * \param[in] t0 Start time [s] of the trace (at the observer).
     * \param[in] dt Duration of one time segment [s] of the observation trace.
     * \param[in] nots Number of time segments.
     * \param[out] ObservationField Electromagnetic field integrals at the observation point.
     */
    void integrateFieldTrace(
	    Vector ObservationPoint,
    	double t0,
    	double dt,
    	int nots,
    	std::vector<ElMagField> *ObservationField);

private:

    /*! Compute the time when a signal emiited from the trajectory point
     * with given index is observed at the ObservationPoint.
     * 
     * @todo watch out! - validity of the index is not checked
     */
    virtual double RetardedTime(int index,
			Vector ObservationPoint);
    
    /*! Compute the field emiited from the trajectory point
     * with given index observed at the ObservationPoint.
     * 
     * @todo watch out! - validity of the index is not checked
     * 
     * This method is private - no methods with direct indexing of
     * the trajectory points should be outside-visible. For the outside
     * world there exists a method with the same name but the
     * observation time as a parameter.
     */
    virtual ElMagField RetardedField(int index,
			     Vector ObservationPoint);

protected:
    
    //! number of trajectory points
    int NP;
    
private:
    
    //! charge in units of ElementaryCharge
    double Charge;
    
    //! mass in unit of the electron rest mass
    double Mass;
    
    //! time in lab frame [s]
    std::vector<double> Time;
    
    //! position in lab frame [m]
    std::vector<Vector> X;
    
    //! dimensionless momentum in lab frame : \f$\beta\gamma = \frac{c p}{mc^2}\f$
    std::vector<Vector> P;
    
    //! acceleration [1/s] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    std::vector<Vector> A;
    
    //! time step for tracking - this will remain constant after being set at the start of tracking
    double dt;
    
    //! charge over mass ratio - this will remain constant during tracking
    double qm;

    //! charge over mass divided by the double time step - this will remain constant during tracking
    double qmt2;
    
    //! \f$u^{i+1}/c\f$ of the Vay algorithm
    Vector VY_p_i1;
    
    //! \f$\gamma^{i+1}/c\f$ of the Vay algorithm
    double VY_gamma_i1;
    
};
