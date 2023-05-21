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

#include "global.h"
#include "fields.h"
#include "fieldtrace.h"
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
 * the laboratory frame. While trajectories are typically computed step-by-step,
 * the class also provides the possibility to pre-allocate the memory for a certain
 * length of the trajectory in order to speed up the access and to reduce the number
 * of memory allocations.
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
    
    /*! Standard constructor:<br>
     * Charge is is given in units of elementary charge, mass im units of electron mass.
     * Memory for a number (nTraj) of trajectory points is pre-allocated.
     */
    ChargedParticle(double charge, double mass, int nTraj);
    
    /*! Copy constructor:<br>
     * Create a replication of a given particle.
     */
    ChargedParticle(const ChargedParticle *origin);
    
    /*! Constructor from buffer:<br>
     *  Take the data stored with serialize() and recreate this
     *  as a new particle including initTrajectory().
     *  The resulting particle will have one trajectory point and no pre-allocated memory.
     *  There is no check of buffer size.
     */
    ChargedParticle(double *buffer);
    
    /*! Constructor from buffer including trajectory:<br>
     *  Take the data stored with serializeTraj() and recreate this
     *  as a new particle including all trajectory data.
     *  There is no check of buffer size.
     *  The number of trajectory points must be known beforehand
     *  as well as the pre-allocated trajectory length.
     */
    ChargedParticle(double *buffer, int nTraj, int nPre);
    
    /*! Destructor: free all memory. */
    virtual ~ChargedParticle();
    
    /*! Remove all trajectory information from the particle.
     */
    void clearTrajectory();

    /*! Pre-allocate a number of trajectory points.
     *  Note - this must be one more than the number of tracking steps
     */
    void preAllocate(int nTraj);

    /*! Return the number of points in the trajectory */
    int getNP() { return NP; };

    /*! Return the number of pre-allocated points in the trajectory */
    int getNpre() { return Npre; };

    /*! Return the charge of the particle in C */
    double getCharge() { return Charge*ElementaryCharge; };
    
    //! Return the time of the particle at the last half-step
    double getTime() { return t_current; };
    
    //! Return the position of the particle at the last half-step
    Vector getPosition() { return X_current; };
    
    //! Return the momentum of the particle at the last half-step
    Vector getMomentum() { return P_current; };
    
    //! Return the acceleration of the particle at the last half-step
    Vector getAccel() { return A_current; };

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
    
    /*! Return the size of a buffer (number of doubles)
     *  that can hold all particle data including all trajectory points.
     *  This buffer size is needed to call serializeTraj()
     */
    int TrajBufSize();
    
    /*! Copy all current information about the particle into one buffer
     *  of double type. This can be used to transfer the particle to a different
     *  MPI node. Only the number of trajectory points is not serialized
     *  and must be known beforehand to restore the particle from the buffer.
     *  The serialized properties include:
     *  @item double Charge
     *  @item double Mass
     *  @item double dt
     *  @item double qm
     *  @item double qmt2
     *  @item double VY_p_i1
     *  @item double VY_gamma_i1
     *  and for all trajectory points:
     *  @item double Time
     *  @item double[3] X
     *  @item double[3] P
     *  @item double[3] A
     *  The buffer must have a size given by TrajBufSize()
     *  There is no check for sufficient buffer size.
     */
    void serializeTraj(double *buffer);

    /*! Return the size of a buffer (number of doubles)
     *  that can hold all current particle data including fields.
     *  This buffer size is needed to call bufferProbe()
     *
     */
    unsigned int getProbeBufferSize() { return 16; };
    
    /*! @brief Buffer particle information.
     * 
     *  Particle coordinates (time,position,momentum,acceleration)
     *  perceived fields (E,B)
     *  1x time, 3x position, 3x momentum, 3x acceleration
     *  3x E, 3x B
     *
     *  The necessary buffer size is given by getProbeBufferSize().
     *  There is no memory check for the buffer size performed.
     * 
     *  @param[in] buffer memory pointer of the buffer
     *  @return the buffer pointer to the next free byte
     */
    double* bufferProbe(double *buffer);

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
    
    /*! Shift the trajectory by a given displacement vector.
     *  The charge is scaled by an additional given factor which is also applied to qm and qmt2.
     */
    void Shift(Vector delta, double f_charge);
    
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

    /*! @brief set current coordinates
     *
     *  This only affects the current particle coordinates.
     *  Trajectory storage is not altered.
     *
     *  Using this method breaks tracking, as additional tracking state variables
     *  are not set properly (VY_p_i1 and VY_gamma_i1).
     *  It is OK to use this to tranfer particle properties from the tracked beam
     *  onto the master beam which is not tracked but used to generate
     *  interaction fields and log files.
     * 
     *  @param[in] time time stamp of the trajectory point
     *  @param[in] pos position
     *  @param[in] mom momentum
     *  @param[in] acc accelearation
     */
    void setCurrentPoint(double time, Vector pos, Vector mom, Vector acc);
    
    /*! @brief add trajectory point
     * 
     *  One trajectory point with the given coordinates is pushed
     *  onto the back of the existing paritcle trajectory.
     * 
     *  @param[in] time time stamp of the trajectory point
     *  @param[in] pos position
     *  @param[in] mom momentum
     *  @param[in] acc accelearation
     */
    void storeTrajectoryPoint(double time, Vector pos, Vector mom, Vector acc);
    
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
     * as defined in the FieldTrace object.
     * 
     * This method sums up the field generated by the particle over all time steps.
     * Fields emitted before the start of the trajectory are not included.
     * The methods RetardedTime() and RetardedField() are called
     * to sequentially process all time steps stored in the trajectory.
     * All field contributions are just added to the field referenced
     * by ObservationField. The caller may choose to start with a zero field
     * or to use this to sum up the fields of multiple particles.
     *
     * This method is supposed to conserve the time integral of the fields.
     * The methods steps over the time-steps of the trajectory
     * and distributes the observed fields into the FieldTrace according
     * to its arrival time. The sampling points of the FieldTrace are
     * considered the centers of dt long buckets. The field of the source
     * is linearly interpolated and distributed onto te buckets accordingly.
     * This way short field spikes falling into the interior of one time step
     * are not overlooked as it would be in a simple sampling algorithm.
     * 
     * \param[in] ObservationPoint The position [m] of the observer.
     * \param[out] trace Electromagnetic field integrals at the observation point.
     */
    void integrateFieldTrace(
	    Vector ObservationPoint,
        FieldTrace *trace);

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
    
    //! the number of trajectory points that can be held in pre-allocated memory
    int Npre;

    //! charge in units of ElementaryCharge
    double Charge;
    
    //! mass in unit of the electron rest mass
    double Mass;
    
    //! current time in lab frame [s]
    //! this is at the half-step position of the Vay tracking algorithm
    double t_current;
    
    //! current position in lab frame [m]
    //! this is at the half-step position of the Vay tracking algorithm
    Vector X_current;
    
    //! current dimensionless momentum in lab frame : \f$\beta\gamma = \frac{c p}{mc^2}\f$
    //! this is at the half-step position of the Vay tracking algorithm
    Vector P_current;
    
    //! current acceleration [1/s] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    //! this is at the half-step position of the Vay tracking algorithm
    Vector A_current;
    
    //! current electric field [V/m] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    //! this is at the half-step position of the Vay tracking algorithm
    Vector E_current;

    //! current magnetic field [T] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    //! this is at the half-step position of the Vay tracking algorithm
    Vector B_current;

    //! trajectory time in lab frame [s]
    std::vector<double> Time;
    
    //! trajectory position in lab frame [m]
    std::vector<Vector> X;
    
    //! trajectory dimensionless momentum in lab frame : \f$\beta\gamma = \frac{c p}{mc^2}\f$
    std::vector<Vector> P;
    
    //! trajectory acceleration [1/s] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    std::vector<Vector> A;
    
    //! time step for tracking - this will remain constant after being set at the start of tracking
    double dt;
    
    //! charge over mass ratio - this will remain constant during tracking
    double qm;

    //! charge over mass divided by the double time step - this will remain constant during tracking
    double qmt2;
    
    //! \f$u^{i+1}/c\f$ of the Vay algorithm (beta*gamma)
    Vector VY_p_i1;
    
    /*! When one tracking step has been computed the particle is half a time step
     *  ahead of its latest stored trajectory point. This is the relativistic
     *  energy of the particle at that position.
     *  \f$\gamma^{i+1}/c\f$ of the Vay algorithm
     */
    double VY_gamma_i1;
    
};
