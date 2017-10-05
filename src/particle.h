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

/*!
 * \class ChargedParticle
 * \brief Particle with electric charge.
 * @author Ulf Lehnert
 * @date 5.10.2017
 * 
 * This class holds the data of a single charged particle with given mass and charge.
 * For tracking this class provides a number of algorithms to perform a single time-step.
 * The trajectory is not stored, only the previous tracking step is preserved.
 * 
 * All charged particles are surrounded by electromagnetic fields and
 * accelerated particle produce radiation. The fields produced by the particle
 * at a given position in laboratory frame can be computed
 * taking into account the full retardation effect.
 * 
 */
class ChargedParticle
{
    
public:
    
    /*! Default constructor:<br>
     * Charge is set to electron charge, mass to electron mass.
     */
    ChargedParticle();
    
    /*! Standard constructor:<br>
     * Charge is is given in units of elementary charge, mass im units of electron mass.
     */
    ChargedParticle(double charge, double mass);
    
    /*! Copy constructor:<br>
     * Create a replication of a given particle.
     */
    ChargedParticle(const ChargedParticle *origin);
    
    /*! Destructor: free all memory. */
    ~ChargedParticle();
    
    /*! Return the charge of the particle */
    int getCharge();
    
    /*! get/set the present time [s] */
    double getTime();
    void setTime(double t);
    
    /*! get/set the particle position [m] */
    Vector getPosition();
    void setPosition(Vector x);
    
    /*! get/set the particle momentum \f$\beta*\gamma\f$ */
    Vector getMomentum();
    void setMomentum(Vector p);
    
    /*! @brief Setup for tracking the particle using the Vay algorithm.
     * 
     * The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     * The initial position, momentum and time of the particle must have been
     * setup before.
     * 
     * @param tstep the length of the time step.
     * @param field the field through which the particle is tracked
     */
    void InitVay(double tstep,
		 GeneralField* field);

    /*! @brief Perform one tracking step using the Vay algorithm.
     * 
     * The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     * @param field the field through which the particle is tracked
     */
    void StepVay(GeneralField* field);
    
    /*! Compute the time when a signal emitted from the particle
     * is observed at the ObservationPoint.
     */
    double PreviousRetardedTime(Vector ObservationPoint);
    double RetardedTime(Vector ObservationPoint);
	
    /*! Compute the field emitted from the particle
     * observed at the ObservationPoint.
     */
    ElMagField PreviousRetardedField(Vector ObservationPoint);
    ElMagField RetardedField(Vector ObservationPoint);

    /*! Compute the time and field generated the particle
     * at an ObservationPoint.
     */
    ElMagObs PreviousObservation(Vector ObservationPoint);
    ElMagObs Observation(Vector ObservationPoint);

private:
    
    //! charge in units of ElementaryCharge
    double Charge;
    
    //! mass in unit of the electron rest mass
    double Mass;
    
    //! time in lab frame [s]
    double PreviousTime;
    double Time;
    
    //! position in lab frame [m]
    Vector PreviousX;
    Vector X;
    
    //! dimensionless momentum in lab frame : \f$\beta\gamma = \frac{c p}{mc^2}\f$
    Vector PreviousP;
    Vector P;
    
    //! acceleration [1/s] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    Vector PreviousA;
    Vector A;
    
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
