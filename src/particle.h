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
 * @date 7.4.2017
 * 
 * This class holds the trajectory data of a single charged particle with given mass and charge.
 * The trajectory comprises time [s], position [m], momentum [] and acceleration [1/s] in
 * the laboratory frame.
 * 
 * For tracking this class provides a number of algorithms to perform a single time-step.
 * 
 * All charged particles are surrounded by electromagnetic fields and
 * accelerated particle produce radiation. The fields produced by the particle
 * at a given time and position in laboratory frame can be computed
 * taking into account the full retardation effect.
 */
class ChargedParticle
{

public:

    /*! Default constructor:<br>
     * Charge is set to electron charge, mass to electron mass.<br>
     * No trajectory is allocated
     */
    ChargedParticle();

    /*! Copy constructor:<br>
     * A given particle is identically reproduced, thereby creating copies of all trajectory data
     */
    ChargedParticle(const ChargedParticle* Part);

    /*! Destructor: free all memory of trajectory storage */
    ~ChargedParticle();

    /*! Return the number of points in the trajectory */
    int GetNP();

    /*! Return the time [s] for one point of the trajectory */
    double TrajTime(int step);

    /*! Return the position [m] for one point of the trajectory */
    Vector TrajPoint(int step);
    
    /*! Return the momentum \f$\beta*\gamma\f$ for one point of the trajectory */
    Vector TrajMomentum(int step);

    // track the particle through a given field
    // use the most simple Euler algorithm
    void TrackEuler(int Nstep,       // number of timesteps
                    double tstep,    // time step size
                    Vector X0,       // initial position
                    Vector P0,       // initial momentum
                    Lattice* field);

    // track the particle through a given field
    // using a leap-frog algorithm
    void TrackLF(int Nstep,       // number of timesteps
                 double tstep,    // time step size
                 Vector X0,       // initial position
                 Vector P0,       // initial momentum
                 Lattice* field);

    // track the particle through a given field
    // using the Vay algorithm
    void TrackVay(int Nstep,       // number of timesteps
                  double tstep,    // time step size
                  Vector X0,       // initial position
                  Vector P0,       // initial momentum
                  Lattice* field);

    // translate a given particle trajectory
    void Translate(Vector R);

    // mirroring a given particle trajectory
    // on a plane y=MirrorY
    // (this also reverses the charge)
    void MirrorY(double MirrorY);

    // electric field radiated by the particle
    // at a given observation point at time t in lab frame
    // retardation is properly accounted for
    ElMagField RetardedField(double time, Vector ObservationPoint);

    /*! Write all information including the trajectory data into an SDDS file.
     * 
     *	 write data on a particle to particle basis
     * 
     *	 use sddsquery trajectory.sdds to view the format of dataset
     *	 
     *	 routine returns values with following meaning:
     *	 
     *	0  ->  Successfully Written the file\n
     *	1  ->  Error in Initializing the Output dataset \n
     *	2  ->  Error in Defining the parameters \n
     *	3  ->  Error in Defining Columns describing the data to be followed \n
     *	4  ->  Error in Writing the layout of the data structure in the sdds file \n
     *	5  ->  Error in Starting a New Page of the file \n
     *	5  ->  Error in setting the values of the parameters \n
     *	6  ->  Error in setting the row values i.e. the data belonging to column \n
     *	7  ->  Error in Writing the page that was successfully initialized \n
     *	8  ->  Error in Terminating the data flow to the sdds file \n
     * 
     */
    int WriteSDDS(const char *filename);
    
private:

    //! number of trajectory points
    int NP;
    
    //! charge in units of ElementaryCharge
    double Charge;
    
    //! mass in unit of the electron rest mass
    double Mass;
    
    //! time in lab frame [s]
    double* Time;
    
    //! position in lab frame [m]
    Vector* X;
    
    //! dimensionless momentum in lab frame : \f$\beta\gamma = \frac{c p}{mc^2}\f$
    Vector* P;
    
    //! acceleration [1/s] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    Vector* A;

};
