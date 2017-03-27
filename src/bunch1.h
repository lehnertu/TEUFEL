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
#include "particle.h"
#include "vector.h"
#include <tuple>
using namespace std;

/*!
    \class Bunch
    \brief Bunching Tracking class
 
    @author Ulf Lehnert, Vipul Joshi
    @date 10.2.2017
    \todo: Read SDDS Input, modify tracking routine
 */

class Bunch {

public:
    /*!

	default constructor: creates single particle resting at origin .
	
	the charge and mass of particles is equal to that of electron 
	
	*/
    Bunch();

    /*
	
	standard constructor:		

	reads a file containing information on
        time,x,y,z,px,py,pz of the particles
        each row of the file must contain above information
	Number of particles equals NP
	must be equal to number of rows
	Charge and mass of the particles have to be provided
	in units of electron's mass and charge

	*/
    Bunch(const char* filename, int NP, int charge, int mass);

    /*!

	copy constructor:

	Create the copy of the bunch

	*/
    Bunch(const Bunch* bunch);

    /*!

	 get the number of particles;

	*/
    int getNOP();

    /*!

	get the number of time steps each particle in bunch undergoes	
	*/
    int getNOTS();

    /*!

	charged particle containing all the particles and its information

	*/
    ChargedParticle* b;

    /*!

	get total charge of the bunch in terms of electron's charge
	
	*/
    int getCharge();

    /*

	get total mass of the bunch of particles in terms of electron's mass

	*/
    int getMass();

    /*!

	Track the bunch through lattice fields using the Euler algorithm 
	interaction fields included
	*/
    void Track_Euler(int NOTS, double tstep, Lattice* field);

    /*!

	Track the bunch through the laatice fields using the Vay Algorithm

	To do: Track individual particles in particle routine
	*/
    void Track_Vay(int NOTS, double tstep, Lattice* field);

    /*!

	returns the radiation field obtained at a time defined by t
	at some observation point defined by the Vector Robs
	for the complete bunch
	*/
    tuple<Vector, Vector> RadiationField(Vector Robs, double t);

    /*!

	 write SDDS file to output file named "trajectory.sdds

	 write data on a particle to particle basis

	 use sddsquery trajectory.sdds to view the format of dataset
	 
	 routine returns values with following meaning:
	 
	0  ->  Successfully Written the file\n
	1  ->  Error in Initializing the Output dataset \n
	2  ->  Error in Defining the parameters \n
	3  ->  Error in Defining Columns describing the data to be followed \n
	4  ->  Error in Writing the layout of the data structure in the sdds file \n
	5  ->  Error in Starting a New Page of the file \n
	5  ->  Error in setting the values of the parameters \n
	6  ->  Error in setting the row values i.e. the data belonging to column \n
	7  ->  Error in Writing the page that was successfully initialized \n
	8  ->  Error in Terminating the data flow to the sdds file \n

	*/
    int WriteSDDSTrajectory();

    /*!

	 write SDDS file to output file named "time-trajectory.sdds

	 writes data for every time step in different page

	 use sddsquery trajectory.sdds to view the format of dataset
	 
	 routine returns values with following meaning:
	 
	0  ->  Successfully Written the file \n
	1  ->  Error in Initializing the Output dataset \n
	2  ->  Error in Defining the parameters \n
	3  ->  Error in Defining Columns describing the data to be followed \n
	4  ->  Error in Writing the layout of the data structure in the sdds file \n
	5  ->  Error in Starting a New Page of the file  \n
	5  ->  Error in setting the values of the parameters \n
	6  ->  Error in setting the row values i.e. the data belonging to column \n
	7  ->  Error in Writing the page that was successfully initialized \n
	8  ->  Error in Terminating the data flow to the sdds file \n

	*/
    int WriteSDDSTime();

private:
    //Number of Particles in the bunch
    int NOP;

    //Charge of the particles in units of electron's charge
    int Charge;

    //Mass of the particles in units of electron mass
    int Mass;

    //filename of the distrubtion file
    const char* file;

    //read the beam profile file and set initial values
    void LoadBeamProfile(const char* filename, const ChargedParticle* part);

    //check for file layout
    int FileCheck(const char* filename, int NP);

    //time step for trajectory integration
    double TIMESTEP;

    //total integration time for the routine
    double TotalTime;

    // Total Charge of the bunch
    double Qtot;

    //Total Mass of the bunch
    double Mtot;

    //total interaction field seen by particle identified by ParticleID
    //fields are calculated with Observation Point==ParticleID
    //time = Laboratory Time or the iteration time
    //stepnumber will always be one less than the iteration number
    tuple<Vector, Vector> MutualField(int stepnumber, int ParticleID, double t);

    //allocate memory to the every particle for stroing the trajectory points
    //arrays will have size NOTS
    void InitializeTrajectory(int NOTS);

    //number of time steps every particle moves;
    int NT;

    // to be used while initializing the bunch
    void setNOP(int NP);
};
