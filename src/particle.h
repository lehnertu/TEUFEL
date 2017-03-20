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

#include "externalfield.h"
#include "vector.h"
#include <tuple>



using namespace std;

/*!
    \class ChargedParticle
    \brief 3D-Particle class
 
    @author Ulf Lehnert, Vipul Joshi
    @date 17.3.2017
    \todo make components private
 */


class ChargedParticle
{

  public:

    /*!
	Default Constructor

	Creates a Single Particle with Charge=-1,Mass=1 and no trajectory points
    */
     ChargedParticle();

    /*!

	Default Constructor

	Create Single Particle with No Trajectory Details
    */ 
    ChargedParticle(int charge, int mass);


    /*!

	Copy Constructor

    */
    ChargedParticle(const ChargedParticle *Part);


    /*!

	Destructor

    */
    ~ChargedParticle();


    /*!

	Set the number of trajectory points the particle can have \n
	Allocate memory for the time,position,momentum,acceleration vectors
    */
    void init( int TrajLength, Vector X0, Vector P0, double T0);    
    

    /*!


	Track the particle through a given field \n
        Use the most simple Euler algorithm 
    */
    void TrackEuler
      (  int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field);

    /*!

	 track the particle through a given field
         Using a leap-frog algorithm
    */
    void TrackLF
      (  int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field);


    /*!

	 Track the particle through a given field
         using the Vay algorithm
    */
    void TrackVay
      (  int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field);

    /*!
	The Vay algorithm defines velocities (momenta) at integer time steps
	positions and field are computed at half-step points.
	We store all quantities at the half-step points.
	The stored velocity is beta_h * gamma_h = u_h / c
	u_h := u^(i+1/2)
	This routine Track the particle through a given "general" field
        using the Vay algorithm for a single step
	it saves the following  variables
	for future trajectories: \n
	1 -> qm (charge to mass ratio)\n
	2 -> t2 (half step timestep value - equal to initial time)\n
	3 -> qmt2(variable used multiple times in Vay Algorithm)\n
	4 -> t_h (half step time value)\n
	5 -> p_h (half step momnetum vector)\n
	6 -> x_h (half step position vector)\n
	6 -> gamma_h (half step relativistic factor) \n
	7 -> beta_h (beta value at half time step)\n
	8 -> dp_dt (acceleration of particle)\n
	9 -> p_i1 ("previous half step" momentum)\n
	10 -> beta_i1 ("previous half step" momentum)\n
	11 -> gamma_i1 ("previous half step" relativistic factor)\n
	 
    */
     void StepVay(
         int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
	 double T0,
         Lattice *field );

    /*! 

	electric field and magnetic field radiated by the particle
        at a given observation point at time t in lab frame.
        retardation is properly accounted for
    */
    tuple<Vector,Vector> RetardedEField(double time, Vector ObservationPoint);


    /*!

	translate a given particle trajectory
    */
    void Translate(Vector R);


    /*!

	 mirroring a given particle trajectory
         on a plane y=MirrorY
        (this also reverses the charge)
    */
    void MirrorY(double MirrorY);

    /*!

	Return the ith trajectory position

    */
    Vector TrajPoint(int i);


    /*!

	Return the time value for ith trajectory point

    */
    double TrajTime(int i);


    /*!

	Return the momentum vector for the ith trajectory point

    */    
    Vector TrajMomentum(int i);
    
    
    /*!


	Return the acceleration vector experienced at ith trajectory point

    */
    Vector TrajAccel(int i);

    /*!

	return the integer value of the Charge
    */
    int getCharge();


    /*!

	get the particle's mass in terms of electron's mass
    */
    int getMass();

 
    /*!

	Return the number of points in the track

    */
     int GetNOTS();
   

    


  private:

    int NOTS;			// number of trajectory points
    int Charge;			// charge in units of ElementaryCharge
    int Mass;			// mass in unit of the electron rest mass
    double *Time;		// time in lab-frame [s]
    Vector *X;			// position in lab frame [m]
    Vector *P;                  // momentum in lab frame : c p = beta gamma mc²
                                // dimensionless in units of mc²
    Vector *A;                  // acceleration in lab frame a = d/dt(p*c)
				// in unit of 1/s (scaled by mc²)

    int counter = 0;		// if counter = 0; then it initializes trajector for the particle
    double qm,t2,qmt2,t_h,gamma_h,gamma_i1;
    Vector x_h,p_h,beta_h, dp_dt,beta_i,p_i1,beta_i1;
};

