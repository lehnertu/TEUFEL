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

#ifndef BUNCH_H
#define BUNCH_H

#include "externalfield.h"
#include "particle.h"
#include "vector.h"
#include <tuple>
using namespace std;

// a class for tracking bunch of charged particles in external and mutual interaction fields

class Bunch
{

  public:
	//default constructor. creates single particle resting at origin .
	// the charge and mass of particles is equal to that of electron 
	Bunch();

	//reads a file containing information on
        //time,x,y,z,px,py,pz of the particles
        //each row of the file must contain above information
	//Number of particles equals NP
	//must be equal to number of rows
	//Charge and mass of the particles have to be provided
	//in units of electron's mass and charge
	Bunch(const char *filename,int NP,int charge,int mass);
	
	//set and get the number of particles;
	int getNOP();
	

	//get the number of time steps each particle in bunch undergoes	
	int getNOTS();



	//charged particle
	ChargedParticle *b;	


	//Track the bunch through lattice fields
	//interaction fields included
	void Track_Euler(int NOTS, double tstep, Lattice *field);
	void Track_Vay(int NOTS, double tstep, Lattice *field);
	


	// returns the radiation field obtained at a time
	// at some observation point
	// for the complete bunch

	tuple<Vector,Vector>RadiationField(Vector Robs, double t);
	
        
    

  private:
	//Number of Particles in the bunch
	int NOP;
	
	//Charge of the particles in units of electron's charge							
	int Charge;

	//Mass of the particles in units of electron mass								
	int Mass;

	//filename of the distrubtion file							
	const char *file;	

	//read the beam profile file and set initial values						
	void LoadBeamProfile(const char *filename,const ChargedParticle *part);

	//check for file layout
	int FileCheck(const char *filename, int NP);

	//total interaction field seen by particle identified by ParticleID
	//fields are calculated with Observation Point==ParticleID
	//time = Laboratory Time or the iteration time
	//stepnumber will always be one less than the iteration number
	tuple<Vector,Vector>MutualField(int stepnumber, int ParticleID, double t);

	//allocate memory to the every particle for stroing the trajectory points
	//arrays will have size NOTS
	void InitializeTrajectory(int NOTS);

	//number of time steps every particle moves;
	int NT;

	// to be used while initializing the bunch
	void setNOP(int NP);

	
   
};

#endif
